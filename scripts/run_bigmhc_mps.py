"""Device-only Apple GPU routing; official BigMHC prediction and weights stay unchanged.

Run with the pinned upstream src on PYTHONPATH. The official CLI accepts CPU/CUDA
only; this wrapper redirects its acceleration hook to MPS without altering the
seven-member ensemble, model forward passes, batching, or output serialization.
"""

import predict
import torch
from bigmhc import BigMHC


def accelerate_mps(model, devices):
    if devices:
        raise ValueError("MPS wrapper requires the official CLI's -d cpu argument")
    return model.to("mps")


if __name__ == "__main__":
    if not torch.backends.mps.is_available():
        raise RuntimeError("Apple MPS is unavailable; use the ordinary CPU runner")
    BigMHC.accelerate = staticmethod(accelerate_mps)
    predict.main()
