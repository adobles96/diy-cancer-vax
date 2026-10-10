// Run with Playwright installed and CHROME_PATH set when using system Chrome.
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const {chromium} = require('playwright');

async function main() {
  const summaries = JSON.parse(fs.readFileSync(process.argv[2] || 'outputs/six-patient-validation.json', 'utf8'));
  const output = path.resolve('outputs/validation');
  fs.mkdirSync(output, {recursive: true});
  const browser = await chromium.launch({headless: true, executablePath: process.env.CHROME_PATH || undefined});
  const checks = [];
  for (const summary of summaries) {
    const page = await browser.newPage({viewport: {width: 1440, height: 1000}});
    const errors = [], network = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.route('**/*', route => {
      const url = route.request().url();
      if (url.startsWith('file:') || url.startsWith('data:')) return route.continue();
      network.push(url);
      return route.abort();
    });
    await page.goto(pathToFileURL(path.join(summary.run_directory, 'report.html')).href, {waitUntil: 'load', timeout: 90000});
    await page.waitForFunction(() => document.querySelector('#historical-chart .main-svg'), {timeout: 30000});
    const originalRows = await page.locator('#epitopes-table tbody tr').count();
    if (originalRows > 40) throw new Error('Large table is not paginated');
    await page.locator('#epitopes-table input').fill('no-such-sequence');
    if (await page.locator('#epitopes-table tbody tr').count()) throw new Error('Epitope filter did not respond');
    await page.locator('#epitopes-table input').fill('');
    await page.locator('#cell-filter').selectOption('CD8');
    await page.locator('#time-filter').selectOption('week 16');
    await page.locator('#result-filter').selectOption('negative');
    const negativeRows = await page.locator('#assays-table tbody tr').count();
    if (!negativeRows) throw new Error('Negative assay rows not shown');
    await page.locator('#result-filter').selectOption('not determined');
    const unknownRows = await page.locator('#assays-table tbody tr').count();
    if (!unknownRows) throw new Error('Unknown assay rows not shown');
    await page.locator('#result-filter').selectOption('');
    if (summary.selection.selected > 1) {
      await page.locator('#sequence-select').selectOption({index: 1});
      const name = await page.locator('#sequence-select option:checked').innerText();
      const detail = await page.locator('#sequence-detail').innerText();
      if (!detail.includes(name.split(' · ')[1])) throw new Error('Sequence selection did not respond');
    }
    await page.locator('#stages button').last().click();
    if (!(await page.locator('#stage-detail').innerText()).includes('after selection')) throw new Error('Pipeline control did not respond');
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.screenshot({path: path.join(output, `html-patient-${summary.patient}.png`)});
    if (summary.patient === 1) {
      await page.locator('#sequences').screenshot({path: path.join(output, 'html-sequence-detail.png')});
      await page.locator('#tumor-table').screenshot({path: path.join(output, 'html-tumor-evidence.png')});
    }
    if (errors.length || network.length) throw new Error(JSON.stringify({errors,network}));
    checks.push({patient: summary.patient, offline: true, javascript_errors: errors, external_requests: network,
                 filter_negative_rows: negativeRows, filter_unknown_rows: unknownRows, initial_table_rows: originalRows,
                 sequence_selection: 'passed', pipeline_buttons: 'passed', graphs: 'rendered'});
    await page.close();
  }
  await browser.close();
  fs.writeFileSync(path.join(output, 'html-checks.json'), JSON.stringify(checks, null, 2));
  console.log(JSON.stringify(checks, null, 2));
}
main().catch(error => {console.error(error); process.exitCode = 1;});
