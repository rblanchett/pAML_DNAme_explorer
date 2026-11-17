from playwright.sync_api import sync_playwright
import json, pathlib
url='http://127.0.0.1:8765/'
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True,executable_path='/usr/bin/chromium',args=['--no-sandbox','--disable-dev-shm-usage'])
 page=browser.new_page(viewport={'width':1600,'height':1000},device_scale_factor=1,accept_downloads=True)
 errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 root=pathlib.Path(__file__).resolve().parent.parent
 html=(root/'index.html').read_text().replace('<script defer src="data/atlas.js"></script><script defer src="assets/app.js"></script>','')
 page.set_content(html,wait_until='load',timeout=80000)
 page.add_style_tag(content=(root/'assets/style.css').read_text())
 page.add_script_tag(content=(root/'data/atlas.js').read_text())
 page.add_script_tag(content=(root/'assets/app.js').read_text())
 page.wait_for_selector('#nodeLayer .node',timeout=30000)
 print('initial title',page.locator('#graphTitle').inner_text())
 assert page.locator('#nodeLayer .node').count()==23
 assert '50,013' in page.locator('#datasetStats').inner_text()
 print('initial edges',page.locator('#edgeLayer .edge').count(),'rows',page.locator('#evidenceBody tr').count())
 page.screenshot(path='/mnt/data/pAML_Methylation_Atlas_Full_v2/tests/screenshot_promoter.png',full_page=True)
 # region changes and gene select
 page.select_option('#region','Gene body')
 print('gene body rows',page.locator('#tableCount').inner_text())
 page.select_option('#region','All')
 page.locator('#geneSearch').fill('WT1')
 page.locator('#suggestions button').filter(has_text='WT1').first.click()
 print('WT1 all regions',page.locator('#graphSubtitle').inner_text())
 assert 'WT1' in page.locator('#graphTitle').inner_text()
 if page.locator('#edgeLayer .edge').count():
  page.locator('#edgeLayer .edge').first.locator('path').last.dispatch_event('click')
  print('clicked arrow',page.locator('#details h3').first.inner_text())
  assert 'Clinical bridge context' in page.locator('#details').inner_text()
 page.locator('#inspectPair').click()
 assert 'Clinical bridge context' in page.locator('#details').inner_text()
 print('clinical pair',page.locator('#details h3').first.inner_text())
 page.locator('#evidenceBody tr').first.click()
 print('site',page.locator('#details h3').first.inner_text())
 assert 'All preserved suite labels' in page.locator('#details').inner_text()
 page.locator('#probeSearch').fill('cg08664529')
 page.locator('#goProbe').click()
 assert 'cg08664529' in page.locator('#graphTitle').inner_text()
 print('exact CpG',page.locator('#graphSubtitle').inner_text())
 with page.expect_download(timeout=30000) as d:
  page.locator('#exportCsv').click()
 dl=d.value
 assert dl.suggested_filename.endswith('.csv')
 print('export',dl.suggested_filename)
 page.locator('details.advanced summary').click()
 # Exercise alternate contrast, island, flags, group filters and 13C public CSV.
 page.select_option('#contrast','IP')
 page.select_option('#contrast','PN')
 page.select_option('#contrast','IN')
 page.select_option('#signature','core')
 page.select_option('#signature','all')
 page.select_option('#island','Open sea')
 page.select_option('#island','All')
 page.select_option('#sourceFilter','nmf6')
 page.select_option('#sourceFilter','All')
 with page.expect_download(timeout=30000) as d2:
  page.locator('#exportClinicalCsv').click()
 assert d2.value.suggested_filename.endswith('Public.csv')
 print('13C export',d2.value.suggested_filename)
 print('errors',errors)
 assert not errors, errors
 page.screenshot(path='/mnt/data/pAML_Methylation_Atlas_Full_v2/tests/screenshot_cpg.png',full_page=True)
 browser.close()
