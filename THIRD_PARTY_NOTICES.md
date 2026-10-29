# Third-party sources

Source workbooks and UCSC protein sequences are not bundled in Git or in the Python package. `src/neovax/resources/sources.json` records upstream download locations and exact hashes. Consult the source publishers' terms before redistributing their material.

NetMHCpan is accessed through the IEDB hosted API. This repository does not contain its source, binaries, weights or a redistribution license. The original historical study used NetMHCpan 2.4; the application explicitly requests binding-affinity version 4.1.

Optional DeepImmuno files are downloaded from the pinned upstream repository into ignored `data/outputs/models/`. The upstream MIT license is downloaded alongside them as `LICENSE` and preserved. Its copyright notice and permission text apply to those files. Model files are not redistributed with this package.

Offline reports bundle Plotly JavaScript obtained through the installed Plotly Python package, retaining the upstream license header. Plotly.js is MIT-licensed. Other Python dependencies retain their own package licenses; see the locked dependency metadata and upstream distributions.

The exact reference-proteome search uses the BSD-3-Clause-licensed `pyahocorasick` package, installed as a dependency. Reference sequences remain separately downloaded, checksum-pinned inputs.
