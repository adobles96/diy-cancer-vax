# DIY Cancer Vaccine

Explore neoantigen selection using published data from six NeoVax patients, or reproduce the charts and tables from the companion blog post.

[Read Alejandro’s Substack](https://substack.com/@alejandro6)

## Start here

Install Git and **Python 3.11**, then:

```sh
git clone https://github.com/adobles96/diy-cancer-vax.git
cd diy-cancer-vax
./run
```

The launcher installs dependencies automatically. A simple terminal menu offers:

- **Run an individual patient** — choose patient 1–6, accept the defaults or change the method and peptide settings, then see the selected long peptides, driving epitopes, scores, genes, variants and historical vaccine matches.
- **Reproduce the blog results** — run all six patients with NetMHCpan, DeepImmuno and their ensemble; save recall@20 and ROC-AUC charts, peptide/tumor recognition tables, and an HTML report.

When source data is missing, you choose where to store it. The program downloads approximately 255 MB, shows its progress and remembers the directory. DeepImmuno’s separate environment and model are prepared automatically when needed. The first run needs internet access and can take a while; subsequent runs reuse downloaded files and predictions.

## Results

```text
data/
├── raw/                 # Published workbooks and reference proteins
└── outputs/
    ├── runs/<run>/      # Patient report, shortlist CSV and one run record
    ├── blog/            # Blog report, charts and tables
    ├── cache/           # Reusable predictions
    └── models/          # Downloaded model files
```

Each patient run contains just `report.html`, `proposed_vaccine.csv` (sequences, annotations, scores and historical matches), and `run.json` (settings, provenance and detailed results). The report lets you explore candidate rankings and assay evidence without separate intermediate files.

The terminal prints the report location when finished. Data and generated outputs are excluded from Git. Patient reports are self-contained; keep the blog report with its charts and tables.

Patient runs default to a 500 nM pMHC binding affinity cutoff. Blog reproduction uses fixed settings: no binding cutoff, top 20 **mutation targets**, and approximate HLA substitutions where DeepImmuno lacks exact support. Metric definitions, model details and limitations are included in the generated reports.

This is an exploratory analysis of published research, not instructions for preparing or administering a vaccine. Sources: [Ott et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5577644/) and [Hu et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8273876/).
