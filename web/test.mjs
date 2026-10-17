import assert from 'node:assert/strict';
import {createServer} from 'node:http';
import {readFile,mkdir} from 'node:fs/promises';
import {resolve,extname} from 'node:path';
import {pathToFileURL} from 'node:url';
import {chromium,devices} from 'playwright';

await mkdir('test-output',{recursive:true});
const server=createServer(async (req,res) => {
  try {
    const file=resolve('site','.'+decodeURIComponent(new URL(req.url,'http://localhost').pathname));
    assert(file.startsWith(resolve('site')));
    const target=file===resolve('site') ? resolve(file,'index.html') : file;
    const mime={'.html':'text/html','.js':'text/javascript','.webmanifest':'application/manifest+json','.png':'image/png','.ico':'image/vnd.microsoft.icon'};
    res.setHeader('Content-Type',mime[extname(target)]||'text/plain');
    res.end(await readFile(target));
  } catch {res.writeHead(404);res.end('Not found');}
});
await new Promise(done=>server.listen(0,'127.0.0.1',done));
const url=`http://127.0.0.1:${server.address().port}/index.html`;
let browser;
try {
  browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,
    args:['--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  const errors=[], external=[];
  const mobile=await browser.newContext({...devices['Pixel 7'], viewport:{width:448,height:998}});
  const page=await mobile.newPage();
  page.on('pageerror',error=>errors.push(error.message));
  page.on('request',request=>{if(!request.url().startsWith('http://127.0.0.1:'))external.push(request.url());});
  await page.goto(url);
  await page.waitForFunction(()=>window.trebViewer && document.getElementById('status').textContent.startsWith('Offline ready'));
  const faviconUrl=await page.locator('link[rel="icon"]').getAttribute('href');
  assert(faviconUrl.startsWith('./favicon.ico?v='),'The page must explicitly request the versioned icon');
  const iconResult=await mobile.request.get(new URL(faviconUrl,url).href);
  assert.equal(iconResult.status(),200);
  assert.equal((await iconResult.body()).subarray(0,4).toString('hex'),'00000100');
  assert(await page.evaluate(async icon => !!(await (await caches.open(window.TREB_CACHE)).match(new URL(icon,location.href))),faviconUrl),
    'The favicon must be cached for offline use');
  assert.equal(await page.locator('#energies .energy-row').count(),7);
  const measurements=await page.evaluate(()=>({range:trebViewer.data.maxRangeFeet,release:trebViewer.data.releaseTime,angle:trebViewer.data.launchAngle,
    viewport:document.getElementById('viewport').getBoundingClientRect().width,overflow:document.documentElement.scrollWidth>innerWidth}));
  assert(Math.abs(measurements.range-2177.325)<.03);
  assert(Math.abs(measurements.release-.683125)<1e-5);
  assert(Math.abs(measurements.angle-30)<1e-8);
  assert.equal(await page.evaluate(()=>trebViewer.data.slingWeightLb),11/16);
  assert.equal(await page.evaluate(()=>trebViewer.data.pumpkinWeightLb),8);
  assert((await page.locator('#details').textContent()).includes('9 ft 6.5 in sling'));
  assert.equal(await page.evaluate(()=>trebViewer.data.armTipHardwareWeightLb),9/16);
  assert(!measurements.overflow);
  await page.screenshot({path:'test-output/phone-cocked.png',fullPage:true});
  // Equal slider movement must mean equal elapsed time across the dense
  // release region, not equal numbers of exported samples.
  assert.equal(Number(await page.locator('#time').getAttribute('max')),2);
  const samples=[];
  const start=Math.floor((await page.evaluate(()=>trebViewer.data.releaseTime))/.01)*.01-.03;
  for(const time of Array.from({length:7},(_,i)=>start+i*.01)) {
    samples.push(await page.evaluate(time=>{
      const slider=document.getElementById('time');slider.value=time;
      slider.dispatchEvent(new Event('input'));
      return {slider:Number(slider.value),time:trebViewer.data.time[trebViewer.index],index:trebViewer.index};
    },time));
  }
  for(let i=1;i<samples.length;i++) {
    assert(Math.abs(samples[i].slider-samples[i-1].slider-.01)<1e-9);
    assert(Math.abs(samples[i].time-samples[i-1].time-.01)<.0021);
  }
  const density=samples.slice(1).map((sample,i)=>sample.index-samples[i].index);
  assert(Math.max(...density)>10*Math.min(...density),
    'The check must exercise the real nonuniform export density');
  await page.locator('#time').focus();
  const arrowStart=Number(await page.locator('#time').inputValue());
  await page.keyboard.press('ArrowRight');await page.keyboard.press('ArrowRight');
  assert(Number(await page.locator('#time').inputValue())>arrowStart,'Fine arrow steps must accumulate between saved samples');
  await page.locator('#release').click();
  assert.equal(await page.locator('#phase').textContent(),'Pumpkin released');
  assert.equal(await page.evaluate(()=>trebViewer.scene.children.find(x=>x.name==='machine').visible),true);
  await page.locator('#frame').click();
  assert.equal(await page.evaluate(()=>trebViewer.scene.children.find(x=>x.name==='machine').visible),false);
  assert.equal(await page.evaluate(()=>trebViewer.scene.children.find(x=>x.name==='tube').visible),true);
  await page.locator('#forces').check();
  await page.locator('#frame').click();
  await page.screenshot({path:'test-output/phone-release.png',fullPage:true});
  const before=await page.evaluate(()=>trebViewer.camera.position.toArray());
  const box=await page.locator('canvas').boundingBox();
  await page.mouse.move(box.x+box.width*.4,box.y+box.height*.5);
  await page.mouse.down(); await page.mouse.move(box.x+box.width*.7,box.y+box.height*.6,{steps:12}); await page.mouse.up();
  await page.waitForTimeout(300);
  assert.notDeepEqual(await page.evaluate(()=>trebViewer.camera.position.toArray()),before);
  await page.locator('#side').click(); await page.locator('#reset').click();
  const stateBefore=await page.evaluate(()=>({index:trebViewer.index,range:document.getElementById('range').textContent,
    energies:document.getElementById('energies').innerHTML,forces:document.getElementById('forces').checked,
    frame:document.getElementById('frame').getAttribute('aria-pressed')}));
  const cameraBefore=await page.evaluate(()=>trebViewer.camera.position.toArray());
  await page.locator('#dimensions').click();
  assert.equal(await page.evaluate(()=>trebViewer.viewMode),'2d');
  assert.equal(await page.evaluate(()=>trebViewer.camera.isOrthographicCamera),true);
  await page.mouse.move(box.x+box.width*.4,box.y+box.height*.5);
  await page.mouse.down(); await page.mouse.move(box.x+box.width*.6,box.y+box.height*.5,{steps:6}); await page.mouse.up();
  await page.waitForTimeout(200);
  const flatDirection=await page.evaluate(()=>{const c=trebViewer.camera,v=new c.position.constructor();c.getWorldDirection(v);return v.toArray();});
  assert(Math.abs(flatDirection[0])<1e-6 && Math.abs(flatDirection[2])<1e-6,'2D drag must pan without rotating');
  assert.deepEqual(await page.evaluate(()=>({index:trebViewer.index,range:document.getElementById('range').textContent,
    energies:document.getElementById('energies').innerHTML,forces:document.getElementById('forces').checked,
    frame:document.getElementById('frame').getAttribute('aria-pressed')})),stateBefore);
  await page.locator('#reset').click();
  await page.screenshot({path:'test-output/phone-2d.png',fullPage:true});
  await page.locator('#dimensions').click();
  assert.equal(await page.evaluate(()=>trebViewer.camera.isPerspectiveCamera),true);
  const cameraRestored=await page.evaluate(()=>trebViewer.camera.position.toArray());
  assert(cameraRestored.every((v,i)=>Math.abs(v-cameraBefore[i])<1e-8),'Returning to 3D must restore its camera');
  // Touch events exercise the same controls Android Chrome uses.
  const session=await mobile.newCDPSession(page);
  const touch={x:box.x+box.width*.4,y:box.y+box.height*.5,id:1};
  const beforeTouch=await page.evaluate(()=>trebViewer.camera.position.toArray());
  await session.send('Input.dispatchTouchEvent',{type:'touchStart',touchPoints:[touch]});
  await session.send('Input.dispatchTouchEvent',{type:'touchMove',touchPoints:[{...touch,x:touch.x+60}]});
  await session.send('Input.dispatchTouchEvent',{type:'touchEnd',touchPoints:[]});
  await page.waitForTimeout(200);
  assert.notDeepEqual(await page.evaluate(()=>trebViewer.camera.position.toArray()),beforeTouch);
  const zoomBefore=await page.evaluate(()=>trebViewer.camera.position.distanceTo(new trebViewer.camera.position.constructor(1.5,0,.35)));
  const second={...touch,id:2,x:touch.x+60};
  await session.send('Input.dispatchTouchEvent',{type:'touchStart',touchPoints:[touch,second]});
  await session.send('Input.dispatchTouchEvent',{type:'touchMove',touchPoints:[{...touch,x:touch.x-25},{...second,x:second.x+25}]});
  await session.send('Input.dispatchTouchEvent',{type:'touchEnd',touchPoints:[]});
  await page.waitForTimeout(200);
  const zoomAfter=await page.evaluate(()=>trebViewer.camera.position.distanceTo(new trebViewer.camera.position.constructor(1.5,0,.35)));
  assert(zoomAfter<zoomBefore,'A pinch-out must zoom in');
  await page.locator('#play').click(); await page.waitForTimeout(250); await page.locator('#play').click();
  assert(await page.evaluate(()=>trebViewer.index>trebViewer.data.releaseIndex),'Playback must advance time');
  await mobile.setOffline(true);
  await page.close();
  const offlinePage=await mobile.newPage();
  offlinePage.on('pageerror',error=>errors.push(error.message));
  await offlinePage.goto(url);
  await offlinePage.waitForFunction(()=>!!window.trebViewer);
  await offlinePage.locator('#release').click();
  assert.equal(await offlinePage.locator('#range').textContent(),'2177.3 ft');
  await offlinePage.goto(url.replace('index.html','treb-2026-offline.html'));
  await offlinePage.waitForFunction(()=>!!window.trebViewer);
  await mobile.close();
  // The download must also work from the filesystem without an HTTP server.
  const local=await browser.newContext({viewport:{width:1200,height:900},offline:true});
  const localPage=await local.newPage();
  const fileRequests=[];
  localPage.on('pageerror',error=>errors.push(error.message));
  localPage.on('request',request=>fileRequests.push(request.url()));
  await localPage.goto(pathToFileURL(resolve('site/treb-2026-offline.html')).href);
  await localPage.waitForFunction(()=>!!window.trebViewer);
  assert((await localPage.locator('link[rel="icon"]').getAttribute('href')).startsWith('data:image/x-icon;base64,'),
    'The download must embed the favicon');
  await localPage.locator('#release').click();
  await localPage.locator('#forces').check();
  await localPage.screenshot({path:'test-output/desktop-release.png',fullPage:true});
  assert.equal(fileRequests.length,1,'The standalone file must not fetch companion assets');
  assert.equal(await localPage.evaluate(()=>trebViewer.renderer.getContext().isContextLost()),false);
  const software=await browser.newContext({viewport:{width:1100,height:900},offline:true});
  await software.addInitScript(()=>{
    const original=HTMLCanvasElement.prototype.getContext;
    HTMLCanvasElement.prototype.getContext=function(kind,...args){
      return kind.startsWith('webgl') ? null : original.call(this,kind,...args);
    };
  });
  const softwarePage=await software.newPage();
  softwarePage.on('pageerror',error=>errors.push(error.message));
  await softwarePage.goto(pathToFileURL(resolve('site/treb-2026-offline.html')).href);
  await softwarePage.waitForFunction(()=>!!window.trebViewer);
  assert.equal(await softwarePage.evaluate(()=>trebViewer.softwareRendering),true);
  assert(await softwarePage.locator('#viewport svg path').count()>100,'Software fallback must draw the machine');
  await softwarePage.locator('#release').click();
  assert.equal(await softwarePage.locator('#range').textContent(),'2177.3 ft');
  await softwarePage.locator('#dimensions').click();
  assert.equal(await softwarePage.evaluate(()=>trebViewer.camera.isOrthographicCamera),true);
  await softwarePage.locator('#dimensions').click();
  await softwarePage.screenshot({path:'test-output/software-release.png',fullPage:true});
  const broken=await browser.newContext();
  await broken.addInitScript(()=>{window.atob=()=>{throw new Error('Injected startup failure');};});
  const brokenPage=await broken.newPage();
  await brokenPage.goto(pathToFileURL(resolve('site/treb-2026-offline.html')).href);
  await brokenPage.waitForFunction(()=>document.getElementById('status').textContent.includes('Injected startup failure'));
  assert.equal(await brokenPage.locator('#details').textContent(),'The viewer could not start.');
  assert.deepEqual(errors,[]); assert.deepEqual(external,[]);
  console.log('PASS: offline file, phone controls, launch metrics, PWA offline cache, software 3D with WebGL disabled, and visible startup errors; no unexpected errors or external requests.');
} finally {
  await browser?.close(); await new Promise(done=>server.close(done));
}
