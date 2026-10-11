"""Explicit NCBI BLAST submission and rate-limited retrieval.

Constructing a client sends no request. submit() and retrieve() are explicit
network operations. The demo only generates request plans.
"""
import json, re, time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from .evidence import validate_dna

ENDPOINT = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def batch_spacers(records, max_bases=1000, max_records=25):
    if max_bases < 1 or max_records < 1:
        raise ValueError("Batch limits must be positive")
    batches, batch, size = [], {}, 0
    for name, sequence in records.items():
        sequence = validate_dna(sequence)
        if not name or any(c.isspace() for c in name) or ">" in name:
            raise ValueError("Invalid FASTA ID")
        if len(sequence) > max_bases:
            raise ValueError("Sequence exceeds upload batch limit")
        if batch and (size + len(sequence) > max_bases or len(batch) >= max_records):
            batches.append(batch)
            batch = {}
            size = 0
        batch[name] = sequence
        size += len(sequence)
    if batch:
        batches.append(batch)
    return batches


def request_parameters(records, email, database="nt", relaxed=False):
    if not records:
        raise ValueError("Cannot submit an empty query")
    if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email):
        raise ValueError("Contact email is required for NCBI requests")
    for name, sequence in records.items():
        if not name or any(c.isspace() for c in name) or ">" in name:
            raise ValueError("Invalid FASTA ID")
        validate_dna(sequence)
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": database,
        "QUERY": "".join(f">{k}\n{validate_dna(v)}\n" for k, v in records.items()),
        "HITLIST_SIZE": "100",
        "email": email,
        "tool": "CRISPRPHIBert",
    }
    if relaxed:
        params.update({"WORD_SIZE": "7", "FILTER": "F", "EXPECT": "1000"})
    return params


class BlastClient:
    def __init__(self, transport=None, clock=time.monotonic):
        self.transport = transport or self._transport
        self.clock = clock
        self.last_request = None
        self.last_rid = {}

    @staticmethod
    def _transport(request):
        with urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8")

    def _send(self, params, rid=None):
        now = self.clock()
        if self.last_request is not None and now - self.last_request < 10:
            raise RuntimeError("Wait at least 10 seconds between NCBI requests")
        if rid is not None and rid in self.last_rid and now - self.last_rid[rid] < 60:
            raise RuntimeError("Wait at least 60 seconds before polling this RID again")
        self.last_request = now
        if rid is not None:
            self.last_rid[rid] = now
        request = Request(
            ENDPOINT,
            data=urlencode(params).encode(),
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": "CRISPRPHIBert/0.1",
            },
            method="POST",
        )
        return self.transport(request)

    def submit(self, records, email, database="nt", relaxed=False):
        response = self._send(request_parameters(records, email, database, relaxed))
        match = re.search(r"\bRID\s*=\s*([A-Za-z0-9-]+)", response)
        if not match:
            raise RuntimeError("BLAST submission returned no request identifier")
        rid = match.group(1)
        self.last_rid[rid] = self.clock()
        return {"rid": rid, "database_requested": database, "query_ids": list(records)}

    def retrieve(self, rid):
        if not re.fullmatch(r"[A-Za-z0-9-]+", rid):
            raise ValueError("Invalid RID")
        response = self._send(
            {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}, rid=rid
        )
        status = re.search(r"Status\s*=\s*(WAITING|FAILED|UNKNOWN)", response)
        if status:
            value = status.group(1)
            if value != "WAITING":
                raise RuntimeError(f"BLAST request status: {value}")
            return {"status": "WAITING", "rid": rid}
        try:
            data = json.loads(response)
        except json.JSONDecodeError as error:
            raise RuntimeError("BLAST response is not JSON2") from error
        if "BlastOutput2" not in data:
            raise RuntimeError("BLAST response has no result records")
        return {"status": "READY", "rid": rid, "result": data}
