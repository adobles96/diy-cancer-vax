# Usage and methodology

## Running the pipeline

`./run` opens a short terminal menu. Choose an individual patient or reproduce the blog results. For a patient, accept the defaults or choose the method, epitope lengths, number of targets, long-peptide length and binding cutoff. DeepImmuno and the ensemble also offer HLA rescue: related-allele substitution when exact model support is missing.

Default patient settings are NetMHCpan 4.1 binding affinity, mutation-containing 9/10-mers, IC50 strictly below 500 nM, and up to 20 targets with one proposed 25-amino-acid peptide per mutation. The best retained short epitope drives each selection. Terminal output includes gene, protein change, long peptide, driving epitope, score and historical mutation membership.

The blog workflow uses **fixed, different settings**: no binding cutoff; minimum NetMHCpan IC50 and maximum DeepImmuno score over retained peptide–HLA pairs per mutation; then within-patient fractional ranks, with averaged ties. The ensemble averages the two mutation ranks with equal weights. It selects 20 distinct genomic mutations per method. HLA rescue is enabled, with substitutions recorded in the report.

```sh
./run --patient 3 --data-root /path/to/data
./run --config configs/patient-1.json
./run --patient 3 --cached-only
./run blog --data-root /path/to/data --output-dir /path/to/results
./run blog --cached-only
./run report data/outputs/runs/<run-directory>
./run --help
```

With no arguments, a terminal is required for the menu. For a noninteractive first run, supply `--data-root`. Saved JSON configurations retain their own paths unless overridden with `--data-root`. Configuration also exposes advanced ranking/aggregation settings; changing these does not reproduce the fixed blog benchmark.

## Setup and storage

The launcher creates `.venv/` and installs the project with standard Python and pip. It does not require uv or environment activation. Python 3.11 supports both models; newer Python works for NetMHCpan-only runs. DeepImmuno automatically gets a separate Python 3.11 / TensorFlow 2.15.1 environment under `<data-root>/outputs/deepimmuno-env`. Python 3.11 must be installed and discoverable as `python3.11` (or already available through the model environment).

Source downloads are verified against pinned SHA-256 hashes. The program reports each file it downloads or reuses. The selected data root is saved in Git-ignored `.neovax.json`. Under it:

- `raw/`: source workbooks and hg19 reference proteins.
- `outputs/cache/`: normalized inputs and model prediction caches.
- `outputs/models/`: pinned model code and weights.
- `outputs/runs/<run>/`: individual patient outputs; each execution has a new directory.
- `outputs/blog/`: main charts, CSV tables, HTML report and detailed cohort provenance. Reproducing again refreshes files in this directory; use `--output-dir` to keep separate versions.

Internet access is needed for missing data, model setup and uncached predictions through the public IEDB API. Requests are serial, batched and cached. Cached-only mode downloads nothing and does not start missing model inference. Missing predictions cause an incomplete result, not silently reduced rankings. Ctrl+C requests cancellation; completed predictions remain cached for the next run. A second Ctrl+C exits immediately.

To reuse source files or legacy predictions:

```sh
./run fetch-data --from-directory /path/to/existing/raw
./run import-cache /path/to/prior/results
```

These explicit maintenance commands use their own `--data-dir` or `--cache-dir` options; they do not read the remembered data root.

Prefer your own environment? Run `python -m pip install -e .`, then `neovax` for the menu or `neovax run --patient 2`. The optional `uv.lock` supports uv users; the standard launcher follows the dependency ranges in `pyproject.toml`.

## Candidate generation and reference filtering

The pipeline independently reconstructs missense substitutions from published variant annotations and pinned hg19 proteins. Silent variants, splice-site variants, indels and unsupported or mismatched annotations are excluded with reasons; their proteins are not guessed. It enumerates mutation-containing short peptides before looking at vaccine membership or assay results.

Each candidate short peptide is checked for an exact full-sequence match anywhere in the pinned reference proteome. Matching peptides are removed before scoring. Other peptides from the same mutation remain eligible; a mutation disappears only if no candidates survive. This does not guarantee that the proposed long peptide lacks normal-protein subsequences. Patient germline proteomes are unavailable.

Reference-filter exports include excluded peptides, matching reference transcripts and zero-based half-open protein coordinates. Historical vaccine rescoring is separate from discovery and cannot reintroduce excluded candidates.

## Metrics and interpretation

**Recall@20** counts actual vaccine components whose mutation appears among the top 20 predicted mutation targets, divided by all actual vaccine components. It measures mutation recovery, not exact recovery of the administered long sequence. The published tables contain 98 patient/IMP components, whereas the paper describes 97 unique neoantigens; reports retain all 98 and disclose the discrepancy. Unreconstructable components remain in this recall denominator.

**Peptide recognition** means any author-called positive week-16 CD4 or CD8 IFN-γ ELISpot peptide-response assay, directly ex vivo or after pre-stimulation. These assays test peptide-pulsed autologous antigen-presenting cells. Minigene, tumor and later-timepoint responses do not enter this label. Recognition tables show positive/included vaccine components and percentages. For model columns, only vaccine components recovered in that method’s top 20 are included. Unvaccinated predictions have no response label.

**ROC-AUC** uses all reconstructable, scored vaccinated mutation targets with usable labels, without a top-20 cutoff: 72 targets, 45 positive and 27 negative in the reproduced cohort. Negatives require explicit negative CD4 and CD8 evidence, not missing measurements. Scores and ensemble ranks are unchanged from discovery. Summary rows give the patient mean and within-patient pair-weighted AUC, which never compares targets from different patients. It is a conditional retrospective benchmark, not accuracy on all tumor mutations.

**Tumor recognition** uses any positive week-16 CD4/CD8 autologous tumor response. Its denominator includes only components with an explicit positive or negative tumor assay: 31 for NeoVax. Missing tests are excluded. Tumor recognition is not a direct tumor-killing measurement.

Both predictors are class-I models; the broad CD4/CD8 endpoint includes biology they do not directly model. DeepImmuno’s score is not a calibrated probability of vaccine response. Numerical HLA rescue is an approximation, not an allele-specific model. Audited gene-suffix repairs require the same protein change, a unique genomic mutation and containment of the historical long peptide in the reconstructed mutant protein.

## Artifacts and reproducibility

Patient runs save configuration, manifest, variant and exclusion tables, candidate and HLA-pair scores, mutation rankings, proposed and actual vaccine components, assay observations, reference-filter evidence, CSV/FASTA shortlists and an HTML report. Source hashes, model identity, status and export hashes are recorded. Patient HTML reports contain their figures; blog HTML uses neighboring chart files and links to CSV/provenance files.

The blog workflow produces `vaccine-recall-at-20` and `peptide-recognition-roc-auc` charts in PNG, SVG and PDF; recognition and ROC-AUC CSV tables; `summary.json`; and detailed `cohort-recall/` and `cohort-auc/` results. It recomputes these from raw data plus validated prediction caches, without needing prior reports or untracked plotting scripts.

The older cohort commands remain available for explicit analyses:

```sh
.venv/bin/python scripts/cohort_recall.py --help
.venv/bin/python scripts/cohort_auc.py --help
```

DeepImmuno uses pinned upstream code and weights. Its published training CSV had no exact candidate-peptide overlap in our audit; this is not proof against every form of related-sequence leakage. NetMHCpan 4.1’s released BA/EL training data had no exact overlap with the 72 evaluated targets; five discovery peptides appeared as EL-negative examples. Audit scripts are retained in `scripts/`. BigMHC is an archived, excluded pilot because its training data overlapped the evaluation targets; neither normal CLI flow uses it.

## Development

```sh
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest
.venv/bin/python -m ruff check src tests scripts
```

Integration tests need local source files and prediction caches. Generated data, environments, models and outputs are ignored by Git.
