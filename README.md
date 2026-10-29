# NeoVax Explorer

A terminal application and offline report for exploring the six published NeoVax patients. It reconstructs missense candidates from exact reference transcripts, scores mutation-containing short peptides, proposes long-peptide boundaries, and then joins the administered vaccine and measured assay evidence.

The application starts from published somatic calls, not raw sequencing. Results are an **exploratory retrospective comparison**, not clinical vaccine recommendations or an independent model benchmark.

## Start

With Python 3.11 or newer installed, three commands run Patient 1 through the pipeline:

```sh
git clone https://github.com/adobles96/diy-cancer-vax.git
cd diy-cancer-vax
./run
```

The launcher detects Python, creates a local `.venv/`, and installs dependencies using standard `venv` and `pip`. No environment activation or package-manager command is required. Dependencies are installed on the first launch and when `pyproject.toml` changes; subsequent launches reuse the environment.

If source data is missing, it asks for a **data directory** (default: `data/`). It creates `raw/` and `outputs/` underneath that directory, downloads approximately 255 MB of pinned workbooks and reference proteins, and runs the pipeline. Existing valid inputs and prediction caches are reused. The directory choice is remembered in Git-ignored `.neovax.json`; `--data-root` changes it. The first prediction run contacts the public IEDB service. Progress, the proposed shortlist, and report paths appear in the terminal.

```sh
./run --patient 3                  # Run another patient
./run tui                          # Open the interactive interface
./run --data-root /path/to/data     # Choose storage without a prompt
./run --cached-only                # Verify local data and use cached predictions only
./run --help                       # List CLI commands
```

For noninteractive first runs, supply `--data-root`. Saved configuration files keep their own paths unless explicitly overridden with `--data-root`. The launcher runs from the repository directory, so its relative paths are consistent. On Windows, use `py -3.11 scripts/bootstrap.py` in place of `./run`.

Prefer your own environment? Install with `python -m pip install -e .` and use `neovax run` or `neovax tui` directly; the same data setup applies. `uv` remains optional: `uv sync --locked` and `uv run neovax run` use the checked-in lockfile. The default pip launcher follows the dependency ranges in `pyproject.toml`; it does not enforce that lockfile. No other workspace is required.

If the workbooks have already been downloaded, copy and verify them without downloading again:

```sh
./run fetch-data --from-directory /path/to/existing/data/raw
```

Local data lives under **`data/`**, with source inputs separate from generated files:

```text
data/
├── raw/          # Pinned source workbooks and reference protein sequences
└── outputs/      # Generated results, caches, models, logs and validation artifacts
    ├── cache/    # Normalized source catalog and prediction caches
    ├── runs/     # Individual patient runs and HTML reports
    ├── cohort-recall/
    ├── cohort-auc/
    └── blog/     # Figures and tables for the blog
```

Both `data/raw/` and `data/outputs/` are ignored by Git. The original workbooks and model weights are fetched separately rather than redistributed in this repository.

## Run without the terminal interface

```sh
./run --config configs/patient-1.json
./run --config configs/patient-6.json --cached-only
./run report data/outputs/runs/<run-directory>
```

There is a saved default configuration for each patient. Paths inside configuration files are relative to the current working directory, or may be absolute. Each execution creates a new run directory; completed run data is not overwritten by subsequent executions. `report` explicitly regenerates the report and figures for an existing run.

Defaults: 9/10-mers, published deduplicated HLA-A/B alleles, pinned `netmhcpan_ba-4.1`, mutant IC50 **strictly below 500 nM**, binding ranking, at most one component per genomic mutation, 25-residue long peptides, and a budget of up to 20. No candidate RNA filter is applied. Setting `threshold` to `null` disables the binding gate.

Optional settings include 8–11-mers for binding, mutation selection with `min`/`mean`/`median`, a per-mutation cap, long lengths 15–30, and a configurable combined-ranking weight. Shortlist exhaustion produces an explicit shortfall rather than ineligible filler.

After rolling-window generation, an always-on **reference-proteome filter** removes any mutant short peptide whose entire amino-acid sequence occurs anywhere in the pinned UCSC hg19 protein reference (`hg19_knownGeneTxPep.txt.gz`, 63,691 transcript protein records). This includes alternate isoforms and other locations in the source protein. Matching is exact and length-specific: a normal 9-mer does not by itself exclude a longer 10-mer. Protein boundaries and ambiguous residues are not crossed or treated as wildcards. This uses the annotated reference proteins, not a six-frame translation of the genome or patient-specific germline data; absent matches do not establish tumor specificity.

Only retained candidates enter discovery scoring, ranking and selection. Their aligned wild-type comparators remain available, and historical vaccine rescores remain separate. Each run exports `reference_filter.json` (counts and source hash) and `reference_filtered_epitopes.csv` (every excluded association with matching transcript IDs and zero-based, half-open protein coordinates). The TUI has a Reference matches tab, and HTML reports show the filter counts and evidence. A mutation is removed from ranking only if no candidate peptides survive.

In the TUI, choose a patient, edit settings, and press **Ctrl+R**. **Ctrl+C** requests cancellation after the current request. The Explore tabs support gene/sequence search, membership/evidence filters, pagination and **Enter** for details. **O** opens the report. **Q** cancels an active run; press it again after exports finish to exit. Narrow terminals can scroll the cards and configuration.

## Caching, requests and incomplete runs

Binding requests go serially to the public [IEDB API](https://tools.iedb.org/main/tools-api/) with the explicit method `netmhcpan_ba-4.1`. Bounded batches, retry limits and rate-limit handling prevent uncontrolled dispatch. Each successful batch is checked for complete sequence-window coverage, coordinates, allele, length and finite scores. Original request FASTA, payload, response and hashes are retained alongside a SQLite lookup index.

Run manifests record `complete`, `partial`, `cancelled` or `failed`. Partial inference cannot produce a complete report. Missing allele scores cannot silently make a peptide eligible. Cancelled and partial runs retain completed cache entries; run the same configuration again to resume. A cached-only run never calls IEDB or starts missing DeepImmuno inference. Changing only ranking or selection settings reuses predictions. CLI partial/cancelled runs return exit code 2 after writing their reports.

Discovery inference and historical-sequence rescoring are separate stages. Historical sequences are scored only after discovery ranking and selection are fixed. Published neoORFs appear in the actual-vaccine browser and comparison-only rescores, never in the reconstructed discovery set. Hosted access to NetMHCpan is not an open-source license; the original paper used version 2.4.

The earlier Patient 1 cache can be imported after validating its original inputs and response hashes:

```sh
./run import-cache /path/to/prior/results
./run --config configs/legacy-patient-1.json
```

That directory must contain `prediction_input.fasta` plus the eight `iedb_*.json` / `.tsv` pairs. Legacy mode has no binding gate and uses mutation selection with `min`: 79 reconstructed candidates, all 11 historical missense targets, and 4/11 target overlap at budget 11. The five absent transcript versions remain visible.

## Optional DeepImmuno

The binding installation does not import TensorFlow. Use a separate Python 3.11 environment:

```sh
python3.11 -m venv data/outputs/deepimmuno-env
data/outputs/deepimmuno-env/bin/python -m pip install \
  tensorflow==2.15.1 'numpy<2' 'pandas<4'
./run setup-deepimmuno --python data/outputs/deepimmuno-env/bin/python
./run --config configs/combined-patient-1.json
```

The wrapper pins upstream commit `df42ac5b6bddfe531268335e2dcb496559cd488b`, source files and weight hashes, and invokes upstream inference in a separate process. It validates every returned pair and blocks upstream's unsupported-allele substitution by default; the cohort script can explicitly enable it as described below. The optional environment was exercised on macOS ARM; other platforms need a compatible TensorFlow installation.

Patient 1 has complete exact allele coverage. The pinned upstream allele table **does not cover all alleles for Patients 2–6**:

| Patient | Unsupported alleles |
| --- | --- |
| 2 | HLA-B*38:01, HLA-B*56:01 |
| 3 | HLA-B*47:01 |
| 4 | HLA-A*25:01, HLA-B*27:02 |
| 5 | HLA-A*66:01 |
| 6 | HLA-A*66:01, HLA-A*01:03 |

The standard pipeline uses binding-only mode for those patients. Unsupported optional configurations preserve available binding data but produce no combined shortlist. The separate cohort analysis below can explicitly enable upstream allele rescue. DeepImmuno supports only 9/10-mers and estimates general CD8 immunogenicity propensity, not a patient's TCR recognition or tumor killing. A released training-CSV exact-sequence overlap audit accompanies optional inference; related sequences, other training inputs, model development and NetMHCpan overlap remain unaudited.

## Cohort recall comparison

After running `./run` and setting up DeepImmuno, cache the six patients' binding predictions and run the cohort scripts. These commands use the default `data/` layout:

```sh
.venv/bin/python scripts/cohort_recall.py --rescue-unknown-hla
.venv/bin/python scripts/cohort_recall.py --rescue-unknown-hla --cached-only
```

To populate the six patients' binding caches first, run `.venv/bin/python scripts/validate_cohort.py --live`.

This creates `data/outputs/cohort-recall/report.html`, plots, per-patient provenance, all mutation rankings, vaccine-component evidence/ranks, and recall at 10–100 targets. The same reference-proteome filter runs before cohort scoring and exports `reference-filter-summary.json` and `reference-filtered-epitopes.csv`. Conditional recall and AUC use the retained, scoreable target set. Outputs remain ignored by Git. The script only reads cached binding predictions; the first command can run missing DeepImmuno inference. Previously cached DeepImmuno batches can supply a filtered subset only when model hashes and actual-to-model allele mappings agree. Without the explicit rescue flag, exact allele support is still required.

The comparison selects **20 distinct mutation targets** from each patient's independently reconstructed missense candidates, representing one vaccine component per mutation. It uses the best IC50 for NetMHCpan and highest DeepImmuno score across mutation-containing 9/10-mers and actual HLA-A/B alleles, with **no binding gate**. The equal-weight ensemble averages the two **mutation-level fractional ranks**, with average ranks for ties. This differs from the standard pipeline's pair-level combination. Final selection ties use genomic mutation ID, without labels.

Recall counts vaccine IMP components recovered by mutation identity. It does not require exact long-peptide sequence recovery. The broad immunogenic endpoint is any positive **week-16 CD4 or CD8 peptide response**, including expanded-cell assays; it does not require minigene, tumor, or killing evidence. Later therapy-associated responses are excluded. Both full-denominator and reconstructable-subset recall are exported. The source tables contain 98 patient/IMP components, versus 97 unique neoantigens in the paper's headline; all 98 table entries are retained, and the unresolved discrepancy is disclosed. Publication gene suffixes can be reconciled only with an exact protein consequence, unique mutation, and full vaccine peptide match to a reconstructed mutant protein.

Rescue follows the pinned upstream numeric-name heuristic, including table-order subtype selection. Every actual-to-model HLA mapping is recorded, and rescued scores have separate cache identities. Binding always uses the actual allele. These results are exploratory: the original vaccine was selected using binding predictions, CD4 responses are outside the direct scope of these class-I models, and model training overlap is incompletely audited.

To measure immunogenicity discrimination within the reconstructable vaccine set, without a top-K cutoff:

```sh
.venv/bin/python scripts/cohort_auc.py
```

This verifies the saved cohort-score file hashes and writes a separate report, ROC plots, metrics, and target-level score/assay audit to ignored `data/outputs/cohort-auc/`. Each genomic mutation counts once. A positive label requires any positive week-16 CD4/CD8 peptide assay. A negative label requires explicit measurements in both cell types and all relevant observations negative; incomplete or ambiguous evidence is excluded. Non-vaccine mutations are not labeled negative. Model and ensemble scores are reused unchanged, including the ensemble's normalization over all discovery candidates.

ROC AUC gives equal scores half credit and is undefined if either class is absent. Summaries include the mean patient AUC and an AUC weighted by the number of within-patient positive–negative pairs; the latter never compares scores between patients. These metrics describe discrimination conditional on prior vaccine selection and reconstruction, with small per-patient sample sizes.

The positive endpoint specifically uses the authors' **IFN-γ ELISpot peptide-response calls**: CD8 after pre-stimulation (Ott supplementary workbook 5 column O and Hu Dataset 4a column O), CD4 ex vivo (Hu Dataset 4b column J), or CD4 after pre-stimulation (column K), all at week 16. We do not recalculate positivity from spot counts. Minigene, tumor, HLA-blocking, tetramer, degranulation, pooled cytokine and later-timepoint results do not enter this label; mutant selectivity and tumor killing are not required. The AUC set and AUC values are independent of the recall budget.

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
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest tests/test_core.py tests/test_network.py
.venv/bin/python scripts/validate_cohort.py --live  # first six-patient inference
.venv/bin/python scripts/validate_cohort.py         # cached-only replay
.venv/bin/python -m pytest --basetemp data/outputs/test-tmp
.venv/bin/python -m ruff check src tests scripts/validate_cohort.py
```

Core tests use synthetic records and mocked HTTP responses, never live IEDB. Integration tests use fetched workbooks and existing caches, including label leakage and six-patient Textual navigation. TUI screenshots and browser-validation results are written to `data/outputs/validation/`. Browser checks use Playwright (`node scripts/validate_html.cjs` after installing Playwright; optionally set `CHROME_PATH`) with all non-file/data requests blocked. They exercise pagination, filters, sequence changes, pipeline controls and plotted figures.

The six-patient run index is `data/outputs/six-patient-validation.json`. Reports are standalone files: move an individual `report.html` anywhere and open it without a server or network connection.

## Optional BigMHC IM comparison

The isolated pilot in `scripts/bigmhc_pilot.py` compares BigMHC **IM** (immunogenicity, not the EL presentation model), DeepImmuno, NetMHCpan, and the two equal-weight rank ensembles on the existing six-patient cohort. It leaves the application defaults and DeepImmuno implementation intact. Every regenerated baseline target, peptide–HLA count, score, cost and rank must match the saved cohort before inference. The pilot retains the reference-proteome filter, mutation-containing 9/10-mers, actual patient HLA-A/B, best-pair mutation aggregation and within-patient ranking. All 19 distinct patient alleles have exact BigMHC pseudosequences; the adapter forbids upstream's fuzzy HLA substitution.

**This is a contaminated retrospective scoring comparison.** The released BigMHC IM training set contains exact retained candidate peptide–HLA pairs for 69 of the 72 evaluated vaccine mutations; training plus validation covers 71/72. Excluding those leaves one positive and no negative, so no uncontaminated AUC can be estimated. Results cannot establish that BigMHC improves independent prediction. The released pair labels also differ from this benchmark's broad week-16 CD4/CD8 peptide-response endpoint. Neither endpoint establishes tumor recognition or killing.

Create a separate Python 3.11 runtime; approximately 2.4 GB of official model files are downloaded into the ignored output directory. The BigMHC Academic License permits noncommercial use, including personal non-profit research; commercial use requires a separate license. Source and all seven ensemble members are pinned to `199e890fa69a63646730bcecea0a520e08ca4f02`, verified against official Git blob hashes, and recorded with SHA-256 checksums. No model files are committed.

```sh
python3.11 -m venv data/outputs/bigmhc-pilot/runtime
data/outputs/bigmhc-pilot/runtime/bin/python -m pip install -r requirements-bigmhc.txt
.venv/bin/python scripts/setup_bigmhc.py
# First save datasets.zip from the released research data as:
# data/outputs/bigmhc-pilot/audit/datasets.zip
.venv/bin/python scripts/audit_bigmhc_overlap.py
.venv/bin/python scripts/bigmhc_pilot.py --device cpu --batch-size 512
```

The saved baseline in `data/outputs/cohort-recall` and its validated binding/DeepImmuno caches must exist first. On Apple Silicon, use `--device mps --batch-size 512` for the optional local GPU runner. It changes only the official acceleration hook; the original model code, weights, seven-sigmoid averaging, batching and output are unchanged. The CLI log still prints `devices: []` because its CPU/CUDA argument parser is reused; `inference/inference.json` records the actual MPS device. CPU and GPU floating-point scores can differ slightly. Run `.venv/bin/python scripts/validate_bigmhc_runtime.py` after an MPS run to compare all mutation-winning pairs plus a reproducible sample across every allele/length stratum against the original CPU CLI, and check top-20/AUC agreement after rescoring the winners. This is a numerical validation subset, not a full CPU replay.

Outputs are isolated under `data/outputs/bigmhc-pilot/`: `report.html`, all target rankings, top-20 membership shifts, per-patient and aggregate recall/recognition/AUC, score distributions, coverage, original inference input/output/log, runtime/model provenance, and `audit/overlap.json`. Recognition fractions apply only to tested vaccine components overlapping a method's top 20; unvaccinated candidates are not treated as negative. AUC uses all 72 reconstructable tested mutations (45 positive, 27 negative), with macro and within-patient pair-weighted summaries. Scores are research outputs, not calibrated response probabilities.
