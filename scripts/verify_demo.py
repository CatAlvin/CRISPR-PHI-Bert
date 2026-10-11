"""Verify reproducible outputs from the labelled synthetic sequence fixture."""

import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
reference = root / "sample-output"
expected = json.loads((reference / "summary.json").read_text(encoding="utf-8"))
actual = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))
if actual != expected:
    raise SystemExit("Demo summary differs from the checked-in synthetic reference")
for name in (
    "expanded_spacers.fasta",
    "flank_fragments.fasta",
    "provenance.json",
    "evidence.tsv",
    "candidate_ranking.json",
    "input_hashes.json",
):
    if (args.output / name).read_bytes() != (reference / name).read_bytes():
        raise SystemExit(f"Demo output differs from reference: {name}")
hashes = json.loads((args.output / "input_hashes.json").read_text(encoding="utf-8"))
for name, expected_hash in hashes.items():
    actual_hash = hashlib.sha256((root / "examples" / name).read_bytes()).hexdigest()
    if actual_hash != expected_hash:
        raise SystemExit(f"Fixture hash differs: {name}")
print(f"Verified synthetic output and {len(hashes)} fixture SHA-256 hashes")
