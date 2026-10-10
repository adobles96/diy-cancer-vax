# NeoVax Explorer

A terminal application and offline report for exploring the six published NeoVax patients. It reconstructs missense candidates from exact reference transcripts, scores mutation-containing short peptides, proposes long-peptide boundaries, and then joins the administered vaccine and measured assay evidence.

The application starts from published somatic calls, not raw sequencing. Results are an **exploratory retrospective comparison**, not clinical vaccine recommendations or an independent model benchmark.

## Start

Python 3.11 or newer is required. With [uv](https://docs.astral.sh/uv/):

```sh
uv sync --locked --extra test
uv run neovax fetch-data
uv run neovax tui
```

Alternatively, create a virtual environment and install `pip install -e '.[test]'`. The source manifest downloads approximately 255 MB of workbooks and reference proteins. Downloads are checked against pinned SHA-256 hashes. No other workspace is required.

If the workbooks have already been downloaded, copy and verify them without downloading again:

```sh
neovax fetch-data --from-directory /path/to/existing/data/raw
```

All data, caches, models, reports, exports, logs and validation artifacts go under **`outputs/`**, which Git ignores. The application never commits them. The original workbooks and model weights are fetched separately rather than redistributed in this repository.

## Run without the terminal interface

```sh
neovax run --config configs/patient-1.json
neovax run --config configs/patient-6.json --cached-only
neovax report outputs/runs/<run-directory>
```

There is a saved default configuration for each patient. Paths inside configuration files are relative to the current working directory, or may be absolute. Each execution creates a new run directory; completed run data is not overwritten by subsequent executions. `report` explicitly regenerates the report and figures for an existing run.

Defaults: 9/10-mers, published deduplicated HLA-A/B alleles, pinned `netmhcpan_ba-4.1`, mutant IC50 **strictly below 500 nM**, binding ranking, at most one component per genomic mutation, 25-residue long peptides, and a budget of up to 20. No candidate RNA filter is applied. Setting `threshold` to `null` disables the binding gate.

Optional settings include 8–11-mers for binding, mutation selection with `min`/`mean`/`median`, a per-mutation cap, long lengths 15–30, and a configurable combined-ranking weight. Shortlist exhaustion produces an explicit shortfall rather than ineligible filler.

In the TUI, choose a patient, edit settings, and press **Ctrl+R**. **Ctrl+C** requests cancellation after the current request. The Explore tabs support gene/sequence search, membership/evidence filters, pagination and **Enter** for details. **O** opens the report. **Q** cancels an active run; press it again after exports finish to exit. Narrow terminals can scroll the cards and configuration.

## Caching, requests and incomplete runs

Binding requests go serially to the public [IEDB API](https://tools.iedb.org/main/tools-api/) with the explicit method `netmhcpan_ba-4.1`. Bounded batches, retry limits and rate-limit handling prevent uncontrolled dispatch. Each successful batch is checked for complete sequence-window coverage, coordinates, allele, length and finite scores. Original request FASTA, payload, response and hashes are retained alongside a SQLite lookup index.

Run manifests record `complete`, `partial`, `cancelled` or `failed`. Partial inference cannot produce a complete report. Missing allele scores cannot silently make a peptide eligible. Cancelled and partial runs retain completed cache entries; run the same configuration again to resume. A cached-only run never calls IEDB or starts missing DeepImmuno inference. Changing only ranking or selection settings reuses predictions. CLI partial/cancelled runs return exit code 2 after writing their reports.

Discovery inference and historical-sequence rescoring are separate stages. Historical sequences are scored only after discovery ranking and selection are fixed. Published neoORFs appear in the actual-vaccine browser and comparison-only rescores, never in the reconstructed discovery set. Hosted access to NetMHCpan is not an open-source license; the original paper used version 2.4.

The earlier Patient 1 cache can be imported after validating its original inputs and response hashes:

```sh
neovax import-cache /path/to/prior/results
neovax run --config configs/legacy-patient-1.json
```

That directory must contain `prediction_input.fasta` plus the eight `iedb_*.json` / `.tsv` pairs. Legacy mode has no binding gate and uses mutation selection with `min`: 79 reconstructed candidates, all 11 historical missense targets, and 4/11 target overlap at budget 11. The five absent transcript versions remain visible.

## Optional DeepImmuno

The binding installation does not import TensorFlow. Use a separate Python 3.11 environment:

```sh
uv venv --python 3.11 outputs/deepimmuno-env
uv pip install --python outputs/deepimmuno-env/bin/python \
  tensorflow==2.15.1 'numpy<2' 'pandas<4'
neovax setup-deepimmuno --python outputs/deepimmuno-env/bin/python
neovax run --config configs/combined-patient-1.json
```

The wrapper pins upstream commit `df42ac5b6bddfe531268335e2dcb496559cd488b`, source files and weight hashes, and invokes upstream inference in a separate process. It validates every returned pair and never permits upstream's unsupported-allele substitution. The optional environment was exercised on macOS ARM; other platforms need a compatible TensorFlow installation.

Patient 1 has complete exact allele coverage. The pinned upstream allele table **does not cover all alleles for Patients 2–6**:

| Patient | Unsupported alleles |
| --- | --- |
| 2 | HLA-B*38:01, HLA-B*56:01 |
| 3 | HLA-B*47:01 |
| 4 | HLA-A*25:01, HLA-B*27:02 |
| 5 | HLA-A*66:01 |
| 6 | HLA-A*66:01, HLA-A*01:03 |

Those patients use binding-only mode. Unsupported optional configurations preserve available binding data but produce no combined shortlist. DeepImmuno supports only 9/10-mers and estimates general CD8 immunogenicity propensity, not a patient's TCR recognition or tumor killing. A released training-CSV exact-sequence overlap audit accompanies optional inference; related sequences, other training inputs, model development and NetMHCpan overlap remain unaudited.

## Scientific boundaries

* Variant identity is patient + published genomic coordinates + reference/alternate alleles. Transcript/position associations remain distinct. Gene alone is never a historical target join.
* Exact transcript versions and reference residues are checked. Ambiguous residues exclude affected windows or long contexts; remote ambiguity does not remove a valid mutation context. Indels and splice variants are explicitly unsupported.
* IC50 (nM), returned percentile rank and immunogenicity are distinct fields. Combined ranking uses fractional ranks over the same eligible unique peptide–HLA population, with averaged ties. Costs are not calibrated probabilities.
* Mutation aggregation operates on distinct eligible short peptides after best-allele reduction. In default epitope selection it is descriptive; in mutation selection it affects ordering. No eligible peptides means no aggregate, not zero.
* Proposed long peptides use the documented centering/clamping rule. No invented residues, padding, unphased combinations, synthesis, dosing, formulation or pooling optimization is included.
* Vaccine membership, selected-target TPM and assay labels join **after** selection. Same-target overlap, exact short containment and exact long-sequence identity are reported separately. Proposed sequences themselves were not necessarily assayed.
* Week-16 CD8 peptide ELISpot, minigene processing, autologous tumor recognition, CD4 assays and later measurements remain separate. Blank, `n.d.`, negative and ambiguous results remain distinct. HLA predictions do not establish experimental restriction. Expanded-cell assays do not inherit an ex-vivo positivity threshold.
* CD8 evaluation uses only tested short peptides found in independently reconstructed, fully scored discovery sequences. Gate failures tie at the bottom; missing/unsupported peptides retain their full tested denominator. Average precision groups ties, and precision at K averages boundary ties. Unselected or untested candidates are never negatives.
* Tumor experiments concern selected reactive cell lines. IFN-γ ELISpot measures secretion rather than direct killing. Post-checkpoint-treatment evidence is not treated as vaccine-only evidence. Figure-only and unextracted supporting observations are labeled explicitly.

## Data and provenance

The package manifest records download URLs and checksums. Assays retain raw cell values, workbook/sheet/cell locations and file hashes; repeated HLA rows do not inflate observation counts. Target joins require exact protein consequence and a unique genomic mutation. Unresolved historical targets remain visible.

| Patient | Published rows | Missense | Reconstructed | Missing exact transcript |
| --- | ---: | ---: | ---: | ---: |
| 1 | 179 | 84 | 79 | 5 |
| 2 | 1,887 | 960 | 899 | 61 |
| 3 | 820 | 395 | 356 | 39 |
| 4 | 1,995 | 1,036 | 960 | 76 |
| 5 | 728 | 350 | 325 | 25 |
| 6 | 3,920 | 1,864 | 1,723 | 141 |

Sources: [Ott et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5577644/), [Hu et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8273876/), and [UCSC hg19 reference proteins](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGeneTxPep.txt.gz). The imports audit all patient 1–6 records in 2017 Table 5 and 2021 Datasets 4a/4b. Table 6 is retained as a source for class-II predictions, which are not CD4 response measurements. Supporting figures are catalogued without claiming exhaustive TCR reanalysis.

## Exports and tests

Each run includes saved configuration, manifest, summary, variant/skip/epitope/pair/mutation tables, reconstructed proteins, proposed and actual vaccine tables, reference scores, assay observations, evaluation JSON, FASTA, and one offline HTML file with embedded Plotly. `figures/` contains SVG and PNG figures. `run_data.json` supports report regeneration.

```sh
uv run pytest tests/test_core.py tests/test_network.py
uv run python scripts/validate_cohort.py --live  # first six-patient inference
uv run python scripts/validate_cohort.py         # cached-only replay
uv run pytest --basetemp outputs/test-tmp
uv run ruff check src tests scripts/validate_cohort.py
```

Core tests use synthetic records and mocked HTTP responses, never live IEDB. Integration tests use fetched workbooks and existing caches, including label leakage and six-patient Textual navigation. TUI screenshots and browser-validation results are written to `outputs/validation/`. Browser checks use Playwright (`node scripts/validate_html.cjs` after installing Playwright; optionally set `CHROME_PATH`) with all non-file/data requests blocked. They exercise pagination, filters, sequence changes, pipeline controls and plotted figures.

The six-patient run index is `outputs/six-patient-validation.json`. Reports are standalone files: move an individual `report.html` anywhere and open it without a server or network connection.
