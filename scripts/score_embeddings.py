"""Exercise the statistical stage using labelled synthetic embeddings."""
import argparse, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from crisprphibert.classifier import ReferenceClassifier

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--input",
    type=Path,
    default=Path(__file__).resolve().parents[1] / "examples/synthetic_embeddings.json",
)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
if args.output.exists():
    parser.error("Output already exists")
data = json.loads(args.input.read_text(encoding="utf-8"))
if data.get("kind") != "synthetic_embeddings":
    parser.error("Expected labelled synthetic embeddings")
model = ReferenceClassifier().fit(data["reference"])
result = {
    "kind": "synthetic_statistical_validation",
    "encoder_executed": False,
    "embedding_dimensions": len(data["reference"][0]),
    "reference_rows": len(data["reference"]),
    "scores": {key: model.score(value) for key, value in data["queries"].items()},
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
