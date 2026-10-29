# Third-party sources

Source workbooks and UCSC protein sequences are not bundled in Git or in the Python package. `src/neovax/resources/sources.json` records upstream download locations and exact hashes. Consult the source publishers' terms before redistributing their material.

NetMHCpan is accessed through the IEDB hosted API. This repository does not contain its source, binaries, weights or a redistribution license. The original historical study used NetMHCpan 2.4; the application explicitly requests binding-affinity version 4.1.

Optional DeepImmuno files are downloaded from the pinned upstream repository into ignored `data/outputs/models/`. The upstream MIT license is downloaded alongside them as `LICENSE` and preserved. Its copyright notice and permission text apply to those files. Model files are not redistributed with this package.

Offline reports bundle Plotly JavaScript obtained through the installed Plotly Python package, retaining the upstream license header. Plotly.js is MIT-licensed. Other Python dependencies retain their own package licenses; see the locked dependency metadata and upstream distributions.

The exact reference-proteome search uses the BSD-3-Clause-licensed `pyahocorasick` package, installed as a dependency. Reference sequences remain separately downloaded, checksum-pinned inputs.

The optional BigMHC pilot downloads unmodified code, pseudosequences and model weights from [KarchinLab/bigmhc](https://github.com/KarchinLab/bigmhc), commit `199e890fa69a63646730bcecea0a520e08ca4f02`, into ignored `data/outputs/bigmhc-pilot/upstream/`. The complete upstream `LICENSE` is preserved there. The BigMHC Academic License permits noncommercial use, including personal non-profit purposes; commercial use requires contacting the licensor. Its redistribution and derivative-work conditions apply to those files. They are not bundled or redistributed with this package. The device-only Apple MPS runner imports the upstream implementation without modifying its files. Citation: Albert et al., *Nature Machine Intelligence* (2023), [doi:10.1038/s42256-023-00694-6](https://doi.org/10.1038/s42256-023-00694-6). The overlap audit uses the authors' [released research data, version 4](https://data.mendeley.com/datasets/dvmz6pkzvb/4), stored separately with source and checksum metadata.
