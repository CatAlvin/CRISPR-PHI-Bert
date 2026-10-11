"""Optional DNABERT-S encoder from a reviewed local model snapshot."""
from pathlib import Path


class LocalDNABERTEncoder:
    def __init__(self, model_directory, allow_custom_code=False):
        import torch
        from transformers import AutoTokenizer, AutoModel

        directory = Path(model_directory).absolute()
        if not directory.is_dir():
            raise ValueError("Model path must be a local directory")
        self.torch = torch
        self.tokenizer = AutoTokenizer.from_pretrained(
            str(directory), local_files_only=True, trust_remote_code=allow_custom_code
        )
        self.model = (
            AutoModel.from_pretrained(
                str(directory),
                local_files_only=True,
                trust_remote_code=allow_custom_code,
            )
            .to("cpu")
            .eval()
        )

    def embed(self, sequence):
        from .evidence import validate_dna

        tokens = self.tokenizer(
            validate_dna(sequence), return_tensors="pt", truncation=True, max_length=512
        )
        with self.torch.inference_mode():
            result = self.model(**tokens)
        hidden = (
            result.last_hidden_state
            if hasattr(result, "last_hidden_state")
            else result[0]
        )
        mask = (
            tokens.get("attention_mask", self.torch.ones_like(tokens["input_ids"]))
            .unsqueeze(-1)
            .to(hidden.dtype)
        )
        return (
            ((hidden * mask).sum(1) / mask.sum(1).clamp(min=1))
            .detach()
            .cpu()
            .numpy()[0]
        )
