import json
import threading
import webbrowser
from dataclasses import replace
from pathlib import Path
from typing import ClassVar

from rich.text import Text
from textual import on, work
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import (
    Button,
    Checkbox,
    Collapsible,
    DataTable,
    Footer,
    Header,
    Input,
    Label,
    RichLog,
    Select,
    Static,
    TabbedContent,
    TabPane,
)

from neovax.config import Config
from neovax.data import load_reference, patient_data
from neovax.iedb import BindingCache, plan_batches
from neovax.peptides import enumerate_epitopes, epitope_sequences, reconstruct
from neovax.pipeline import run_pipeline
from neovax.reference_filter import filter_reference_peptides


class DetailScreen(ModalScreen):
    BINDINGS: ClassVar = [("escape", "dismiss", "Close")]

    def __init__(self, record, evidence):
        super().__init__()
        self.record = record
        self.evidence = evidence

    def compose(self) -> ComposeResult:
        with VerticalScroll(id="detail"):
            yield Label(f"{self.record.get('gene', 'Record')} · sequence and evidence")
            sequence = (
                self.record.get("sequence")
                or self.record.get("peptide")
                or self.record.get("immunizing_sequence")
            )
            if sequence:
                if self.record.get("reference_sequence"):
                    yield Static(f"Reference: {self.record['reference_sequence']}", markup=False)
                text = Text(sequence)
                offset = self.record.get("mutation_offset")
                if isinstance(offset, int) and 0 <= offset < len(sequence):
                    text.stylize("bold white on dark_red", offset, offset + 1)
                yield Static(text)
            if self.record.get("driving_epitope"):
                yield Static(
                    f"Driving epitope: {self.record['driving_epitope']} · {self.record['best_allele']}\nIC50: {self.record['ic50_nm']} nM · ranking cost: {self.record['ranking_cost']}",
                    markup=False,
                )
            yield Label(
                "Applicable observations (target response does not prove this proposed sequence was tested)"
            )
            yield Static(
                "\n".join(
                    f"{r['cell_type']} · {r['assay_type']} · {r['timepoint']} · {r['pre_stimulation']}: {r['result']} ({r['sequence']})"
                    for r in self.evidence
                )
                or "No joined observations.",
                markup=False,
            )
            with Collapsible(title="Technical provenance and original source cells"):
                yield Static(
                    json.dumps({"record": self.record, "observations": self.evidence}, indent=2), markup=False
                )
            yield Button("Close", id="close-detail")

    @on(Button.Pressed, "#close-detail")
    def close_detail(self):
        self.dismiss()


class NeoVaxApp(App):
    TITLE = "NeoVax Explorer"
    SUB_TITLE = "Published candidates · predictions · evidence"
    CSS = """
    Screen { background: #10272f; }
    TabPane { padding: 1 2; }
    .patient { width: 1fr; min-width: 22; min-height: 13; height: auto; padding: 1; margin-bottom: 1; border: solid #366a6d; }
    #patient-cards { height: 1fr; }
    .patient-row { height: auto; }
    .patient Button { width: 100%; }
    #configure-scroll { height: 1fr; }
    Input, Select { margin-bottom: 1; }
    Label { margin-top: 1; width: 1fr; height: auto; }
    #configure-scroll > Static { width: 1fr; height: auto; }
    #run, #open-report { margin-top: 1; }
    #progress-log { height: 1fr; border: solid #366a6d; }
    #results-tabs { height: 1fr; }
    DataTable { height: 1fr; }
    #filters { height: auto; }
    #filter-text { width: 2fr; }
    #membership-filter { width: 1fr; }
    #result-summary { height: auto; margin-bottom: 1; }
    #paging { height: 3; }
    #page-label { width: 1fr; content-align: center middle; }
    #detail { width: 90%; height: 90%; margin: 2 4; padding: 1 2; background: #183b43; border: thick #46978f; }
    """
    BINDINGS: ClassVar = [
        ("q", "quit", "Quit"),
        ("ctrl+r", "start_run", "Run"),
        ("ctrl+c", "cancel_run", "Cancel"),
        ("o", "open_report", "Open report"),
    ]

    def __init__(self, config=None):
        super().__init__()
        self.config = config or Config()
        self.result = None
        self.running_pipeline = False
        self.cancel_event = threading.Event()
        self.card_data = {}
        self.page = 0
        self.visible_records = {}
        current_lengths = ",".join(str(x) for x in self.config.lengths)
        self.length_options = [("9 and 10 residues", "9,10"), ("8–11 residues (binding only)", "8,9,10,11")]
        if current_lengths not in {value for _, value in self.length_options}:
            self.length_options.append((f"Saved lengths: {current_lengths}", current_lengths))

    def compose(self) -> ComposeResult:
        yield Header()
        with TabbedContent(initial="patients", id="main-tabs"):
            with TabPane("1 · Patient", id="patients"):
                yield Static(
                    "Choose a published patient. HLA-C is not available. Coverage includes exact transcript matches only."
                )
                with VerticalScroll(id="patient-cards"):
                    for pair in [(1, 2), (3, 4), (5, 6)]:
                        with Horizontal(classes="patient-row"):
                            for patient in pair:
                                with Vertical(classes="patient"):
                                    yield Button(f"Patient {patient}", id=f"patient-{patient}")
                                    yield Static("Loading source coverage…", id=f"card-{patient}")
            with TabPane("2 · Configure", id="configure"):
                with VerticalScroll(id="configure-scroll"):
                    yield Static(f"Patient {self.config.patient}", id="chosen-patient")
                    yield Label("Short peptide lengths")
                    yield Select(
                        self.length_options,
                        value=",".join(str(x) for x in self.config.lengths),
                        allow_blank=False,
                        id="lengths",
                    )
                    yield Label("Binding eligibility: mutant IC50 below this nM (blank disables)")
                    yield Input(
                        "" if self.config.threshold is None else str(self.config.threshold), id="threshold"
                    )
                    yield Label("Ranking method")
                    yield Select(
                        [(x.capitalize(), x) for x in ["binding", "immunogenicity", "combined"]],
                        value=self.config.ranking,
                        allow_blank=False,
                        id="ranking",
                    )
                    yield Static(
                        "Binding: lower nM is better. DeepImmuno: higher propensity is better. Combined: weighted fractional ranks, not a calibrated probability."
                    )
                    yield Label("Vaccine budget (maximum components)")
                    yield Input(str(self.config.budget), id="budget", type="integer")
                    yield Label("Long peptide length (15–30 residues)")
                    yield Input(str(self.config.long_length), id="long-length", type="integer")
                    yield Checkbox(
                        "Cached predictions only (no network)",
                        value=self.config.cached_only,
                        id="cached-only",
                    )
                    with Collapsible(title="Advanced selection and model settings", collapsed=True):
                        yield Label("Selection unit")
                        yield Select(
                            [("Epitope", "epitope"), ("Mutation", "mutation")],
                            value=self.config.selection_unit,
                            allow_blank=False,
                            id="selection-unit",
                        )
                        yield Label("Mutation aggregation (changes selection only in mutation mode)")
                        yield Select(
                            [(x, x) for x in ["min", "mean", "median"]],
                            value=self.config.aggregation,
                            allow_blank=False,
                            id="aggregation",
                        )
                        yield Label("Combined weight on binding (0–1)")
                        yield Input(str(self.config.alpha), id="alpha")
                        yield Label("Maximum components per mutation")
                        yield Input(str(self.config.per_mutation_cap), id="mutation-cap", type="integer")
                        yield Label("Pinned DeepImmuno manifest path (optional isolated model)")
                        yield Input(self.config.deepimmuno_manifest or "", id="deep-manifest")
                    yield Button("Run pipeline", variant="primary", id="run")
                    yield Static("", id="config-error", markup=False)
            with TabPane("3 · Progress", id="progress"):
                yield Static(
                    "Completed batches survive cancellation. A later run with the same inputs resumes from cache."
                )
                yield RichLog(id="progress-log", wrap=True, markup=False)
                yield Button("Cancel after current request", id="cancel")
            with TabPane("4 · Explore", id="explore"):
                yield Static("Run the pipeline to explore results.", id="result-summary", markup=False)
                with Horizontal(id="filters"):
                    yield Input(placeholder="Filter gene, sequence or evidence", id="filter-text")
                    yield Select(
                        [
                            ("All rows", "all"),
                            ("Selected", "selected"),
                            ("Historical target", "historical"),
                            ("Positive assay", "positive"),
                            ("Negative assay", "negative"),
                            ("Missing evidence", "missing"),
                        ],
                        value="all",
                        allow_blank=False,
                        id="membership-filter",
                    )
                with TabbedContent(id="results-tabs"):
                    for title, name in [
                        ("Long peptides", "proposed"),
                        ("Short epitopes", "epitopes"),
                        ("Mutations", "mutations"),
                        ("Actual vaccine", "actual"),
                        ("Assays", "assays"),
                        ("Reference matches", "reference_filtered"),
                    ]:
                        with TabPane(title, id=f"tab-{name}"):
                            yield DataTable(id=f"table-{name}", cursor_type="row")
                with Horizontal(id="paging"):
                    yield Button("Previous", id="previous")
                    yield Static("", id="page-label")
                    yield Button("Next", id="next")
                yield Button("Open offline HTML report", id="open-report")
        yield Footer()

    def on_mount(self):
        self.load_cards()

    @work(thread=True)
    def load_cards(self):
        try:
            for patient in range(1, 7):
                variants, alleles, published, _ = patient_data(self.config.data_dir, patient)
                reference = load_reference(self.config.data_dir)
                proteins, _skipped = reconstruct(variants, reference)
                epitopes, _filtered, filter_summary = filter_reference_peptides(
                    enumerate_epitopes(proteins, self.config.lengths), reference
                )
                cache = BindingCache(self.config.cache_dir)
                try:
                    batches = plan_batches(
                        epitope_sequences(epitopes),
                        alleles,
                        self.config.lengths,
                        cache.lookup(alleles),
                        self.config.batch_size,
                    )
                finally:
                    cache.close()
                card = f"{', '.join(alleles)}\n{len(proteins)}/{len(variants)} variant rows reconstructed\n{len({r['imp_id'] for r in published})} actual vaccine components\n{filter_summary['unique_peptides_removed']} reference-matching peptides excluded\n{len(batches)} binding requests uncached"
                self.card_data[patient] = card
                self.call_from_thread(self.query_one(f"#card-{patient}", Static).update, card)
        except (OSError, ValueError, KeyError) as error:
            for patient in range(1, 7):
                self.call_from_thread(
                    self.query_one(f"#card-{patient}", Static).update, f"Data unavailable: {error}"
                )

    @on(Button.Pressed)
    def button_pressed(self, event):
        button_id = event.button.id or ""
        if button_id.startswith("patient-"):
            if self.running_pipeline:
                return
            patient = int(button_id.split("-")[1])
            self.config = replace(self.config, patient=patient)
            self.query_one("#chosen-patient", Static).update(
                f"Patient {patient}\n{self.card_data.get(patient, '')}"
            )
            self.query_one("#main-tabs", TabbedContent).active = "configure"
        elif button_id == "run":
            self.action_start_run()
        elif button_id == "cancel":
            self.action_cancel_run()
        elif button_id == "open-report":
            self.action_open_report()
        elif button_id in {"previous", "next"}:
            self.page = max(0, self.page + (1 if button_id == "next" else -1))
            self.populate_results()

    def read_config(self):
        def value(selector):
            return self.query_one(selector, Input).value.strip()

        return replace(
            self.config,
            lengths=tuple(int(x) for x in self.query_one("#lengths", Select).value.split(",")),
            threshold=float(value("#threshold")) if value("#threshold") else None,
            ranking=self.query_one("#ranking", Select).value,
            budget=int(value("#budget")),
            long_length=int(value("#long-length")),
            cached_only=self.query_one("#cached-only", Checkbox).value,
            selection_unit=self.query_one("#selection-unit", Select).value,
            aggregation=self.query_one("#aggregation", Select).value,
            alpha=float(value("#alpha")),
            per_mutation_cap=int(value("#mutation-cap")),
            deepimmuno_manifest=value("#deep-manifest") or None,
        )

    def action_start_run(self):
        if self.running_pipeline:
            return
        try:
            self.config = self.read_config()
        except ValueError as error:
            self.query_one("#config-error", Static).update(str(error))
            self.query_one("#main-tabs", TabbedContent).active = "configure"
            self.notify(str(error), severity="error", timeout=8)
            return
        self.query_one("#config-error", Static).update("")
        self.cancel_event.clear()
        self.running_pipeline = True
        self.query_one("#run", Button).disabled = True
        self.query_one("#main-tabs", TabbedContent).active = "progress"
        self.execute_pipeline()

    @work(thread=True, exclusive=True, group="pipeline")
    def execute_pipeline(self):
        try:
            result = run_pipeline(
                self.config,
                progress=lambda message: self.call_from_thread(
                    self.query_one("#progress-log", RichLog).write, message
                ),
                cancel=self.cancel_event.is_set,
            )
        except (OSError, ValueError, KeyError) as error:
            self.call_from_thread(self.pipeline_failed, str(error))
        else:
            self.call_from_thread(self.pipeline_finished, result)

    def pipeline_failed(self, message):
        self.running_pipeline = False
        self.query_one("#run", Button).disabled = False
        self.query_one("#progress-log", RichLog).write(f"Failed: {message}")

    def pipeline_finished(self, result):
        self.result = result
        self.running_pipeline = False
        self.query_one("#run", Button).disabled = False
        summary = result["summary"]
        self.query_one("#result-summary", Static).update(
            f"Patient {summary['patient']} · {summary['status'].upper()} · {len(result['proposed'])} components · {summary['overlap']['target_overlap']} historical targets\n{summary['run_directory']}\nEnter a row for sequence and evidence details."
        )
        self.page = 0
        self.populate_results()
        self.query_one("#main-tabs", TabbedContent).active = "explore"

    @on(Input.Changed, "#filter-text")
    @on(Select.Changed, "#membership-filter")
    def filters_changed(self):
        self.page = 0
        self.populate_results()

    @on(TabbedContent.TabActivated, "#results-tabs")
    def result_tab_changed(self):
        self.page = 0
        self.populate_results()

    def populate_results(self):
        if not self.result:
            return
        active = self.query_one("#results-tabs", TabbedContent).active.removeprefix("tab-")
        if active not in self.result:
            return
        query = self.query_one("#filter-text", Input).value.lower()
        mode = self.query_one("#membership-filter", Select).value
        selected = {r["mutation_id"] for r in self.result["proposed"]}

        def matches(row):
            if query and query not in json.dumps(row).lower():
                return False
            if mode == "selected":
                return row.get("mutation_id") in selected
            if mode == "historical":
                return bool(row.get("historical_target") or active == "actual")
            if mode in {"positive", "negative"}:
                return row.get("result") == mode
            if mode == "missing":
                return (
                    row.get("result") in {"not reported", "not determined", "ambiguous"}
                    or row.get("all_alleles_scored") is False
                )
            return True

        rows = [r for r in self.result[active] if matches(r)]
        self.page = min(self.page, max(0, (len(rows) - 1) // 200))
        self.visible_records[active] = rows[self.page * 200 : (self.page + 1) * 200]
        columns = {
            "proposed": ["gene", "sequence", "driving_epitope", "ranking_cost", "historical_target"],
            "epitopes": ["gene", "peptide", "best_allele", "ic50_nm", "ranking_cost", "eligible"],
            "mutations": ["gene", "protein_change", "ranking_cost", "eligible_epitopes", "total_epitopes"],
            "actual": ["imp_id", "gene", "immunizing_sequence", "reconstructable"],
            "assays": ["gene", "cell_type", "assay_type", "timepoint", "result", "sequence"],
            "reference_filtered": ["gene", "protein_change", "peptide", "reference_matches"],
        }[active]
        table = self.query_one(f"#table-{active}", DataTable)
        table.clear(columns=True)
        table.add_columns(*columns)
        for index, row in enumerate(self.visible_records[active]):
            table.add_row(
                *(str(row.get(c)) if row.get(c) is not None else "—" for c in columns), key=str(index)
            )
        self.query_one("#page-label", Static).update(
            f"{len(rows)} rows · page {self.page + 1}/{max(1, (len(rows) + 199) // 200)}"
        )
        self.query_one("#previous", Button).disabled = self.page == 0
        self.query_one("#next", Button).disabled = (self.page + 1) * 200 >= len(rows)

    @on(DataTable.RowSelected)
    def row_selected(self, event):
        name = event.data_table.id.removeprefix("table-")
        record = self.visible_records[name][int(event.row_key.value)]
        evidence = [
            r
            for r in self.result["assays"]
            if (record.get("imp_id") and r["imp_id"] == record["imp_id"])
            or (record.get("gene") == r["gene"] and record.get("protein_change") == r["protein_change"])
        ]
        self.push_screen(DetailScreen(record, evidence))

    def action_cancel_run(self):
        self.cancel_event.set()
        if self.running_pipeline:
            self.query_one("#progress-log", RichLog).write(
                "Cancellation requested; preserving completed predictions and exports."
            )

    def action_open_report(self):
        if self.result:
            webbrowser.open((Path(self.result["summary"]["run_directory"]) / "report.html").as_uri())

    def action_quit(self):
        if self.running_pipeline:
            self.action_cancel_run()
            self.query_one("#main-tabs", TabbedContent).active = "progress"
        else:
            self.exit()
