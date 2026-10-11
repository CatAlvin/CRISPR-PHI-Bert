from pathlib import Path
from dataclasses import replace
from tempfile import TemporaryDirectory
import json, unittest
from unittest.mock import patch
from crisprphibert.evidence import *
from crisprphibert.__main__ import run_demo
from crisprphibert.tools import local_blast

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        temp_root = Path(__file__).absolute().parent / ".tmp"
        temp_root.mkdir(exist_ok=True)
        self.temp = TemporaryDirectory(dir=temp_root)
        self.root = Path(self.temp.name)
        if self.root.absolute().parent != temp_root:
            raise ValueError("Test directory escaped the workspace")
        self.hit = BlastHit("h", "q", "s", 100, 6, 0, 0, 1, 6, 3, 8, 1e-6, 40)

    def tearDown(self):
        self.temp.cleanup()

    def file(self, text, name="input.txt"):
        p = self.root / name
        p.write_text(text, encoding="utf-8")
        return p

    def test_fasta_multiline_and_lowercase(self):
        self.assertEqual(
            read_fasta(self.file(">a description\nacgt\nNN\n\n")), {"a": "ACGTNN"}
        )

    def test_fasta_duplicate_ids_rejected(self):
        with self.assertRaises(ValueError):
            read_fasta(self.file(">a\nAC\n>a\nGT\n"))

    def test_fasta_empty_header_rejected(self):
        with self.assertRaises(ValueError):
            read_fasta(self.file(">\nAC\n"))

    def test_fasta_empty_sequence_rejected(self):
        with self.assertRaises(ValueError):
            read_fasta(self.file(">a\n>b\nAC\n"))

    def test_fasta_missing_header_rejected(self):
        with self.assertRaises(ValueError):
            read_fasta(self.file("ACGT\n"))

    def test_invalid_nucleotides_rejected(self):
        with self.assertRaises(ValueError):
            validate_dna("ACTZ")

    def test_reverse_complement_iupac_roundtrip(self):
        sequence = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(reverse_complement(sequence)), sequence)

    def test_repeat_identity_with_insertion(self):
        self.assertAlmostEqual(repeat_identity("ACGT", "ACGGT"), 4 / 5)

    def test_protein_coverage_is_percent(self):
        path = self.file("q\ts\t95\t20\t1\t0\t1\t20\t1\t20\t1e-10\t40\t50\n")
        self.assertEqual(select_homologs(path), {})

    def test_protein_reference_level(self):
        path = self.file("q\ts\t99\t20\t0\t0\t1\t20\t1\t20\t1e-40\t40\t99\n")
        self.assertEqual(select_homologs(path), {"s": "LEVEL-2"})

    def test_expansion_deduplicates_sequences(self):
        db = {
            "p": {
                "array_id": "a",
                "repeat": "ACGT",
                "spacers": ["AAAA", "CCCC", "CCCC"],
            }
        }
        expanded, provenance = expand_spacers(
            {"original": "AAAA"}, {"p": "LEVEL-3"}, db, "ACGT"
        )
        self.assertEqual(expanded, {"original": "AAAA", "a_2": "CCCC"})
        self.assertEqual(provenance["a_2"]["source"], "p")

    def test_expansion_reverse_repeat_orients_spacers(self):
        db = {"p": {"array_id": "a", "repeat": "AGTC", "spacers": ["AACG"]}}
        expanded, _ = expand_spacers({}, {"p": "LEVEL-3"}, db, "GACT", 1)
        self.assertEqual(expanded["a_1"], "CGTT")

    def test_expansion_rejects_repeat_mismatch(self):
        db = {"p": {"array_id": "a", "repeat": "AAAA", "spacers": ["AACG"]}}
        self.assertEqual(expand_spacers({}, {"p": "LEVEL-3"}, db, "CCCC", 1)[0], {})

    def test_blast_wrong_columns_rejected(self):
        with self.assertRaises(ValueError):
            read_blast(self.file("q\ts\t100"), {"q": 6})

    def test_blast_unknown_query_rejected(self):
        with self.assertRaises(ValueError):
            read_blast(EXAMPLES / "spacer_blast.tsv", {})

    def test_blast_duplicate_lines_removed(self):
        line = "q\ts\t100\t6\t0\t0\t1\t6\t3\t8\t1e-6\t40\n"
        self.assertEqual(len(read_blast(self.file(line + line), {"q": 6})), 1)

    def test_blast_invalid_query_coordinates_rejected(self):
        line = "q\ts\t100\t6\t0\t0\t1\t7\t3\t9\t1e-6\t40\n"
        with self.assertRaises(ValueError):
            read_blast(self.file(line), {"q": 6})

    def test_full_query_strict_hit(self):
        self.assertTrue(is_significant(self.hit, 6))

    def test_exactly_one_terminal_base_missing(self):
        self.assertTrue(is_significant(replace(self.hit, qstart=2, mismatches=2), 6))

    def test_two_terminal_bases_missing_not_strict(self):
        self.assertFalse(is_significant(replace(self.hit, qstart=2, qend=5), 6))

    def test_terminal_exception_keeps_evalue_constraint(self):
        self.assertFalse(is_significant(replace(self.hit, qstart=2, evalue=2), 6))

    def test_three_mismatches_fail_terminal_exception(self):
        self.assertFalse(is_significant(replace(self.hit, qstart=2, mismatches=3), 6))

    def test_flank_one_based_and_clipped(self):
        self.assertEqual(flank_sequence(self.hit, "AACCGGTTAACC", 2), "AACCGGTTAA")

    def test_flank_reverse_complement(self):
        self.assertEqual(
            flank_sequence(replace(self.hit, sstart=8, send=3), "AACCGGTTAACC", 2),
            "TTAACCGGTT",
        )

    def test_flank_boundary_without_negative_slice(self):
        self.assertEqual(
            flank_sequence(replace(self.hit, sstart=1, send=2), "ACGTT", 400), "ACGTT"
        )

    def test_flank_out_of_bounds_rejected(self):
        with self.assertRaises(ValueError):
            flank_sequence(self.hit, "ACGT", 2)

    def test_strict_hit_retained_without_scores(self):
        self.assertTrue(merge_evidence([self.hit], {"q": 6}, {})[0]["accepted"])

    def test_score_join_is_alignment_specific(self):
        a = replace(self.hit, hit_id="a", qstart=2, qend=5)
        b = replace(a, hit_id="b", sstart=10, send=13)
        rows = merge_evidence([a, b], {"q": 6}, {"a": 0.95, "b": 0.2})
        self.assertEqual([r["accepted"] for r in rows], [True, False])

    def test_score_exact_threshold_accepted(self):
        hit = replace(self.hit, qstart=2, qend=5)
        self.assertTrue(merge_evidence([hit], {"q": 6}, {"h": 0.9})[0]["accepted"])

    def test_unknown_score_id_rejected(self):
        with self.assertRaises(ValueError):
            merge_evidence([self.hit], {"q": 6}, {"other": 0.9})

    def test_nonfinite_score_rejected(self):
        with self.assertRaises(ValueError):
            merge_evidence([self.hit], {"q": 6}, {"h": float("nan")})

    def test_empty_demo_output_rejected_on_second_run(self):
        out = self.root / "result"
        run_demo(EXAMPLES, out)
        with self.assertRaises(FileExistsError):
            run_demo(EXAMPLES, out)

    def test_demo_dataflow_and_explicit_provenance(self):
        out = self.root / "result"
        summary = run_demo(EXAMPLES, out)
        self.assertEqual(
            (summary["input_spacers"], summary["expanded_spacers"]), (3, 7)
        )
        self.assertEqual(
            (
                summary["strict_alignments"],
                summary["embedding_supported"],
                summary["unresolved"],
            ),
            (3, 1, 2),
        )
        self.assertFalse(summary["encoder_executed"])
        self.assertEqual(len(read_fasta(out / "flank_fragments.fasta")), 6)

    def test_blast_wrapper_checks_exit_and_uses_argument_array(self):
        source = self.file(">q\nACGTTT\n", "query with spaces.fasta")
        output = self.root / "output with spaces.tsv"

        def created(command, **kwargs):
            output.write_text("", encoding="utf-8")

        with patch(
            "crisprphibert.tools.subprocess.run", side_effect=created
        ) as process:
            local_blast(source, self.root / "db with spaces", output, relaxed=True)
        command = process.call_args.args[0]
        self.assertIn(str(source), command)
        self.assertTrue(process.call_args.kwargs["check"])
        self.assertFalse(process.call_args.kwargs["shell"])


if __name__ == "__main__":
    unittest.main()
