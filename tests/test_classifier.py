from pathlib import Path
import importlib.util, json, unittest

AVAILABLE = (
    importlib.util.find_spec("numpy") is not None
    and importlib.util.find_spec("sklearn") is not None
)


@unittest.skipUnless(AVAILABLE, "Optional statistics dependencies are not installed")
class ClassifierTests(unittest.TestCase):
    def setUp(self):
        from crisprphibert.classifier import ReferenceClassifier

        self.cls = ReferenceClassifier
        self.fixture = json.loads(
            (
                Path(__file__).resolve().parents[1]
                / "examples/synthetic_embeddings.json"
            ).read_text()
        )

    def test_near_embedding_scores_above_far_embedding(self):
        model = self.cls().fit(self.fixture["reference"])
        self.assertGreater(
            model.score(self.fixture["queries"]["near_reference"])[
                "compatibility_score"
            ],
            model.score(self.fixture["queries"]["far_reference"])[
                "compatibility_score"
            ],
        )

    def test_insufficient_reference_rejected(self):
        with self.assertRaises(ValueError):
            self.cls().fit([[0, 0], [1, 1]])

    def test_singular_reference_covariance_regularized(self):
        model = self.cls().fit([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]])
        score = model.score([2, 2])["compatibility_score"]
        self.assertTrue(0 <= score <= 1)

    def test_nan_embedding_rejected(self):
        with self.assertRaises(ValueError):
            self.cls().fit([[0, 0], [1, 1], [float("nan"), 2]])

    def test_feature_dimension_mismatch_rejected(self):
        model = self.cls().fit(self.fixture["reference"])
        with self.assertRaises(ValueError):
            model.score([1, 2])

    def test_scoring_before_fit_rejected(self):
        with self.assertRaises(ValueError):
            self.cls().score([1, 2])
