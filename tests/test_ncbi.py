import json, unittest
from urllib.parse import parse_qs
from crisprphibert.ncbi import batch_spacers, request_parameters, BlastClient


class SubmissionTests(unittest.TestCase):
    def test_batch_limits_preserve_all_query_ids(self):
        rows = {f"q{i}": "A" * 30 for i in range(40)}
        batches = batch_spacers(rows, max_bases=300, max_records=8)
        self.assertEqual(sum(len(b) for b in batches), 40)
        self.assertTrue(
            all(len(b) <= 8 and sum(map(len, b.values())) <= 300 for b in batches)
        )

    def test_oversized_query_rejected(self):
        with self.assertRaises(ValueError):
            batch_spacers({"q": "A" * 1001})

    def test_upload_parameters_encode_fasta(self):
        params = request_parameters({"q": "acgt"}, "research@example.org", relaxed=True)
        self.assertEqual(params["QUERY"], ">q\nACGT\n")
        self.assertEqual(params["WORD_SIZE"], "7")

    def test_empty_upload_rejected(self):
        with self.assertRaises(ValueError):
            request_parameters({}, "research@example.org")

    def test_submission_tracks_rid_and_request_body(self):
        calls = []
        client = BlastClient(
            transport=lambda req: (calls.append(req) or "RID = TEST-RID\n"),
            clock=lambda: 0,
        )
        result = client.submit({"q": "ACGT"}, "research@example.org")
        self.assertEqual(result["rid"], "TEST-RID")
        self.assertEqual(parse_qs(calls[0].data.decode())["QUERY"], [">q\nACGT\n"])

    def test_missing_rid_is_explicit_failure(self):
        client = BlastClient(transport=lambda req: "error", clock=lambda: 0)
        with self.assertRaises(RuntimeError):
            client.submit({"q": "ACGT"}, "research@example.org")

    def test_retrieval_pacing_and_ready_result(self):
        now = [0]
        responses = iter(["RID = ABC123", json.dumps({"BlastOutput2": []})])
        client = BlastClient(
            transport=lambda req: next(responses), clock=lambda: now[0]
        )
        client.submit({"q": "ACGT"}, "research@example.org")
        now[0] = 10
        with self.assertRaises(RuntimeError):
            client.retrieve("ABC123")
        now[0] = 60
        self.assertEqual(client.retrieve("ABC123")["status"], "READY")

    def test_waiting_status_preserved(self):
        client = BlastClient(transport=lambda req: "Status=WAITING", clock=lambda: 0)
        self.assertEqual(client.retrieve("ABC123")["status"], "WAITING")
