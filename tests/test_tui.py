from pathlib import Path

import pytest
from textual.widgets import DataTable, Input, TabbedContent

from neovax.config import Config
from neovax.tui import DetailScreen, NeoVaxApp

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not Path("data/outputs/cache/binding.sqlite").exists(), reason="Patient prediction cache required"
    ),
]


@pytest.mark.parametrize("patient", range(1, 7))
async def test_patient_config_run_results_navigation(patient, tmp_path):
    app = NeoVaxApp(Config(patient=patient, cached_only=True, output_dir=str(tmp_path)))
    async with app.run_test(size=(110, 40)) as pilot:
        for _ in range(100):
            await pilot.pause(0.1)
            if len(app.card_data) == 6:
                break
        button = app.query_one(f"#patient-{patient}")
        button.scroll_visible(immediate=True)
        await pilot.pause(0.1)
        await pilot.click(f"#patient-{patient}")
        assert app.query_one("#main-tabs", TabbedContent).active == "configure"
        assert app.config.patient == patient
        app.query_one("#budget", Input).value = "5"
        await pilot.press("ctrl+r")
        assert app.query_one("#main-tabs", TabbedContent).active in {"progress", "explore"}
        for _ in range(180):
            await pilot.pause(0.2)
            if not app.running_pipeline:
                break
        assert app.result is not None
        assert app.result["summary"]["status"] == "complete"
        assert app.query_one("#main-tabs", TabbedContent).active == "explore"
        table = app.query_one("#table-proposed", DataTable)
        assert table.row_count == 5
        table.focus()
        await pilot.press("enter")
        assert isinstance(app.screen, DetailScreen)
        await pilot.press("escape")
        app.query_one("#results-tabs", TabbedContent).active = "tab-reference_filtered"
        await pilot.pause(0.1)
        assert app.query_one("#table-reference_filtered", DataTable).row_count == min(
            200, len(app.result["reference_filtered"])
        )
        assert app.result["manifest"]["reference_filter"]["reference_protein_records"] > 0
        app.query_one("#results-tabs", TabbedContent).active = "tab-assays"
        await pilot.pause(0.1)
        assert app.query_one("#table-assays", DataTable).row_count > 0
        app.query_one("#filter-text", Input).value = "no-such-gene"
        await pilot.pause(0.1)
        assert app.query_one("#table-assays", DataTable).row_count == 0
        app.query_one("#filter-text", Input).value = ""
        await pilot.pause(0.1)
        Path("data/outputs/validation").mkdir(parents=True, exist_ok=True)
        app.save_screenshot(f"tui-patient-{patient}.svg", "data/outputs/validation")


async def test_invalid_configuration_and_narrow_layout(tmp_path):
    app = NeoVaxApp(Config(cached_only=True, output_dir=str(tmp_path)))
    async with app.run_test(size=(64, 28)) as pilot:
        await pilot.pause(0.2)
        await pilot.click("#patient-1")
        app.query_one("#long-length", Input).value = "10"
        await pilot.press("ctrl+r")
        assert not app.running_pipeline
        assert "15–30" in str(app.query_one("#config-error").render())
        app.save_screenshot("tui-narrow.svg", "data/outputs/validation")
