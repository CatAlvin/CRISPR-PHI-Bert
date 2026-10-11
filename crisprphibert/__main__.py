from pathlib import Path
from collections import Counter, defaultdict
from dataclasses import asdict
import argparse, csv, hashlib, html, json
from .evidence import (
    read_fasta,
    write_fasta,
    select_homologs,
    expand_spacers,
    read_blast,
    flank_sequence,
    merge_evidence,
)


def dump(path, value):
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def run_demo(data, outdir, flank=20, threshold=0.9):
    data, outdir = Path(data).absolute(), Path(outdir).absolute()
    if not data.is_dir():
        raise FileNotFoundError(data)
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    config = json.loads((data / "config.json").read_text(encoding="utf-8"))
    if config.get("data_kind") != "synthetic":
        raise ValueError("This demo requires its labelled synthetic fixture")
    original = read_fasta(data / "input_spacers.fasta")
    database = json.loads((data / "homologs.json").read_text(encoding="utf-8"))
    homologs = select_homologs(data / "protein_blast.tsv")
    expanded, provenance = expand_spacers(
        original, homologs, database, config["repeat"]
    )
    lengths = {name: len(seq) for name, seq in expanded.items()}
    hits = read_blast(data / "spacer_blast.tsv", lengths)
    genomes = read_fasta(data / "subjects.fasta")
    fragments = {h.hit_id: flank_sequence(h, genomes[h.subject], flank) for h in hits}
    score_fixture = json.loads(
        (data / "illustrative_scores.json").read_text(encoding="utf-8")
    )
    if score_fixture.get("kind") != "illustrative_manual_scores":
        raise ValueError("Unexpected score provenance")
    evidence = merge_evidence(hits, lengths, score_fixture["scores"], threshold)
    counts = Counter(row["evidence"] for row in evidence)
    subject_support = defaultdict(set)
    for row in evidence:
        if row["accepted"]:
            subject_support[row["subject_id"]].add(row["query_id"])
    ranking = sorted(
        (
            {"subject_id": subject, "unique_spacers": len(ids)}
            for subject, ids in subject_support.items()
        ),
        key=lambda x: (-x["unique_spacers"], x["subject_id"]),
    )
    summary = {
        "data_kind": "synthetic",
        "score_kind": "illustrative_manual_scores",
        "encoder_executed": False,
        "input_spacers": len(original),
        "expanded_spacers": len(expanded),
        "protein_candidates": len(database),
        "protein_hits_selected": len(homologs),
        "alignment_records": len(hits),
        "strict_alignments": counts["strict_alignment"],
        "embedding_supported": counts["embedding_supported"],
        "unresolved": counts["unresolved"],
        "accepted_candidates": len(ranking),
        "flank_each_side_bp": flank,
        "score_threshold": threshold,
    }
    outdir.mkdir(parents=True)
    write_fasta(outdir / "expanded_spacers.fasta", expanded)
    write_fasta(outdir / "flank_fragments.fasta", fragments)
    dump(outdir / "provenance.json", provenance)
    dump(outdir / "summary.json", summary)
    dump(outdir / "candidate_ranking.json", ranking)
    hashes = {}
    for path in sorted(data.iterdir()):
        if path.is_file():
            with path.open("rb") as handle:
                hashes[path.name] = hashlib.file_digest(handle, "sha256").hexdigest()
    dump(outdir / "input_hashes.json", hashes)
    with (outdir / "evidence.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(evidence[0])
            if evidence
            else [
                "hit_id",
                "query_id",
                "subject_id",
                "evidence",
                "compatibility_score",
                "accepted",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(evidence)
    rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(x))}</td>"
            for x in [
                r["hit_id"],
                r["query_id"],
                r["subject_id"],
                r["evidence"],
                r["compatibility_score"]
                if r["compatibility_score"] is not None
                else "—",
            ]
        )
        + "</tr>"
        for r in evidence
    )
    (outdir / "report.html").write_text(
        f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CRISPRPHIBert evidence demo</title><style>body{{font:16px/1.65 system-ui,sans-serif;color:#20322c;max-width:1080px;margin:60px auto;padding:0 24px;background:#f8faf7}}h1{{font-size:38px}}.note{{background:#e7eee3;padding:20px;border-radius:8px}}table{{border-collapse:collapse;width:100%;margin:28px 0}}th,td{{text-align:left;border-bottom:1px solid #ccd5cd;padding:12px 8px}}code{{font-size:14px}}@media(max-width:700px){{table{{font-size:12px}}th,td{{padding:8px 3px}}}}</style><p>CRISPRPHIBert / Sequence evidence</p><h1>From spacer expansion to candidate evidence</h1><p class="note">Synthetic sequences and cached BLAST-format alignments. Scores are illustrative inputs. DNABERT inference was not run. This report verifies data handling and decision rules.</p><p>{len(original)} input spacers → {len(expanded)} unique spacers · {counts['strict_alignment']} strict alignments · {counts['embedding_supported']} score-supported alignments · {counts['unresolved']} unresolved alignments</p><table><thead><tr><th>Alignment</th><th>Spacer</th><th>Candidate</th><th>Evidence</th><th>Score</th></tr></thead><tbody>{rows}</tbody></table><p>Each fragment retains an alignment ID through extraction, scoring and merge. The ranking counts unique spacer support per candidate.</p><p>合成序列与缓存比对用于验证数据流程。示例分数不代表模型结果或宿主感染概率。</p></html>""",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="CRISPRPHIBert reproducible sequence evidence demo"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        required=True,
        help="Run the labelled offline fixture",
    )
    parser.add_argument(
        "--data", type=Path, default=Path(__file__).resolve().parents[1] / "examples"
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=20)
    parser.add_argument("--threshold", type=float, default=0.9)
    args = parser.parse_args()
    if args.flank < 0:
        parser.error("--flank must be non-negative")
    try:
        result = run_demo(args.data, args.outdir, args.flank, args.threshold)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
