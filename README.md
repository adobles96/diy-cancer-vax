# DIY Cancer Vaccine

Explore neoantigen selection using published data from six NeoVax patients, or reproduce the charts and tables from the companion blog post.

[Read Alejandro’s Substack](https://substack.com/@alejandro6) · [Usage and methodology](docs/guide.md)

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
    ├── runs/<run>/      # Patient report, shortlist CSV/FASTA and all scores
    ├── blog/            # Blog report, charts and tables
    ├── cache/           # Reusable predictions
    └── models/          # Downloaded model files
```

The terminal prints the report location when finished. Data and generated outputs are excluded from Git. Patient reports are self-contained; keep the blog report with its charts and tables.

For repeatable, unattended runs:

```sh
./run --patient 2 --data-root ./data
./run blog --data-root ./data
./run blog --cached-only             # No downloads or new predictions
```

Patient defaults use NetMHCpan with a 500 nM binding cutoff. Blog reproduction uses the fixed analysis settings: no binding cutoff, top 20 **mutation targets**, and HLA rescue for unsupported DeepImmuno alleles. See the [guide](docs/guide.md) for definitions and limitations. Python 3.12+ also works for NetMHCpan-only runs; DeepImmuno requires Python 3.11 alongside it.

This is an exploratory analysis of published research, not instructions for preparing or administering a vaccine. Sources: [Ott et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5577644/) and [Hu et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8273876/).
