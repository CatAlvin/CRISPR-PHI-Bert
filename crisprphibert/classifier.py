"""Reference distribution from embeddings, without encoder fine-tuning.

The original PCA/GMM/Mahalanobis formula is retained. The output is an
uncalibrated compatibility score; it is not a probability of infection.
"""
import numpy as np
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture


class ReferenceClassifier:
    def fit(self, embeddings):
        x = np.asarray(embeddings, dtype=float)
        if x.ndim != 2 or x.shape[0] < 3 or x.shape[1] < 2 or not np.isfinite(x).all():
            raise ValueError(
                "Reference needs at least 3 finite embeddings with 2 features"
            )
        self.features = x.shape[1]
        self.pca = PCA(n_components=2, svd_solver="full")
        reduced = self.pca.fit_transform(x)
        self.gmm = GaussianMixture(
            n_components=3,
            covariance_type="tied",
            reg_covar=1e-6,
            random_state=0,
            n_init=1,
        ).fit(reduced)
        self.mean = reduced.mean(axis=0)
        self.inverse_cov = np.linalg.pinv(
            np.cov(reduced, rowvar=False) + np.eye(2) * 1e-6
        )
        return self

    def score(self, embedding):
        if not hasattr(self, "pca"):
            raise ValueError("Fit a reference distribution first")
        x = np.asarray(embedding, dtype=float)
        if x.shape != (self.features,) or not np.isfinite(x).all():
            raise ValueError("Invalid query embedding")
        reduced = self.pca.transform(x[None, :])[0]
        delta = reduced - self.mean
        distance = float(np.sqrt(max(0, float(delta @ self.inverse_cov @ delta))))
        log_density = float(self.gmm.score_samples(reduced[None, :])[0])
        z = 0.7 * log_density - 0.3 * distance
        value = (
            float(1 / (1 + np.exp(-z)))
            if z >= 0
            else float(np.exp(z) / (1 + np.exp(z)))
        )
        return {
            "compatibility_score": value,
            "mahalanobis_distance": distance,
            "gmm_log_density": log_density,
        }
