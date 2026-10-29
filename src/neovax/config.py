import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Config:
    patient: int = 1
    lengths: tuple[int, ...] = (9, 10)
    threshold: float | None = 500.0
    ranking: str = "binding"
    alpha: float = 0.5
    selection_unit: str = "epitope"
    aggregation: str = "min"
    long_length: int = 25
    budget: int = 20
    per_mutation_cap: int = 1
    cached_only: bool = False
    batch_size: int = 128
    data_dir: str = "data/raw"
    cache_dir: str = "data/outputs/cache"
    output_dir: str = "data/outputs/runs"
    deepimmuno_manifest: str | None = None
    rescue_unknown_hla: bool = False

    def __post_init__(self):
        object.__setattr__(self, "lengths", tuple(sorted(set(self.lengths))))
        if self.patient not in range(1, 7):
            raise ValueError("Patient must be 1–6")
        if not self.lengths or not set(self.lengths) <= {8, 9, 10, 11}:
            raise ValueError("Lengths must be a nonempty subset of 8–11")
        if self.threshold is not None and (not math.isfinite(self.threshold) or self.threshold <= 0):
            raise ValueError("Binding threshold must be positive and finite, or null to disable")
        if self.ranking not in {"binding", "immunogenicity", "combined"}:
            raise ValueError("Unknown ranking method")
        if self.selection_unit not in {"epitope", "mutation"} or self.aggregation not in {
            "min",
            "mean",
            "median",
        }:
            raise ValueError("Unknown selection unit or aggregation")
        if not 0 <= self.alpha <= 1 or not 15 <= self.long_length <= 30:
            raise ValueError("Alpha must be 0–1 and long length 15–30")
        if self.budget < 1 or self.per_mutation_cap < 1 or not 1 <= self.batch_size <= 256:
            raise ValueError("Budget and cap must be positive; batch size must be 1–256")
        if self.ranking != "binding" and (not set(self.lengths) <= {9, 10} or not self.deepimmuno_manifest):
            raise ValueError(
                "DeepImmuno ranking requires 9/10-mers and a pinned model manifest; use binding otherwise"
            )

    def to_dict(self):
        return asdict(self)

    @classmethod
    def load(cls, path):
        return cls(**json.loads(Path(path).read_text()))
