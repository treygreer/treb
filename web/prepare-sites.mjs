// Prepare and visually check the separate homepage and simulator destinations.
import {chromium,devices} from 'playwright';
import {createServer} from 'node:http';
import {readFile,mkdir,writeFile} from 'node:fs/promises';
import {resolve,extname} from 'node:path';
import {pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const home=resolve('../.site-deploy/firstinfright');
await mkdir('test-output',{recursive:true});
const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,
 args:['--use-angle=swiftshader','--enable-unsafe-swiftshader']});
let serveLegacy=false;
const server=createServer(async(req,res)=>{
 try {
  const pathname=new URL(req.url,'http://localhost').pathname;
  if(serveLegacy && (pathname==='/' || pathname==='/index.html')) {
    res.setHeader('Content-Type','text/html');
    return res.end('<title>Previous viewer</title><script>navigator.serviceWorker.register("./service-worker.js")</script>');
  }
  if(serveLegacy && pathname==='/service-worker.js') {
    res.setHeader('Content-Type','text/javascript');
    return res.end('self.addEventListener("install",event=>event.waitUntil(caches.open("treb-viewer-migration-test").then(cache=>cache.addAll(["./","./index.html"])).then(()=>self.skipWaiting())));self.addEventListener("activate",event=>event.waitUntil(self.clients.claim()));self.addEventListener("fetch",event=>event.respondWith(caches.match(event.request).then(cached=>cached||fetch(event.request))));');
  }
  const root=pathname.startsWith('/treb/') ? resolve('site') : home;
  const relative=pathname.startsWith('/treb/') ? pathname.slice(5) : pathname;
  const file=resolve(root,'.'+(relative.endsWith('/') ? relative+'index.html' : relative));
  assert(file.startsWith(root));
  res.setHeader('Content-Type',({'.html':'text/html','.js':'text/javascript','.png':'image/png','.webmanifest':'application/manifest+json','.ico':'image/x-icon'})[extname(file)]||'text/plain');
  res.end(await readFile(file));
 }catch{res.writeHead(404);res.end('Not found');}
});
await new Promise(done=>server.listen(0,'127.0.0.1',done));
const base=`http://127.0.0.1:${server.address().port}`;
try {
 const page=await browser.newPage({viewport:{width:1140,height:1120}});
 await page.goto(pathToFileURL(resolve('site/treb-2026-offline.html')).href);
 await page.waitForFunction(()=>!!window.trebViewer);
 await page.locator('#release').click();
 await page.waitForTimeout(100);
 await page.locator('#viewport').screenshot({path:resolve(home,'machine.png')});
 await page.goto(base+'/');
 assert.equal(await page.title(),'First in Fright — Carrboro, NC');
 assert.equal(await page.locator('a.button').getAttribute('href'),'https://treygreer.github.io/treb/');
 await page.screenshot({path:'test-output/home-desktop.png',fullPage:true});
 const mobile=await browser.newContext({...devices['Pixel 7'],viewport:{width:448,height:998}});
 const phone=await mobile.newPage(); await phone.goto(base+'/');
 assert.equal(await phone.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false);
 await phone.screenshot({path:'test-output/home-phone.png',fullPage:true});
 // Verify relative assets, scope and offline reopening at the /treb/ project path.
 await phone.goto(base+'/treb/');
 await phone.waitForFunction(()=>document.getElementById('status').textContent.startsWith('Offline ready'));
 assert.equal(await phone.evaluate(()=>new URL(document.querySelector('link[rel="manifest"]').href).pathname),'/treb/manifest.webmanifest');
 await mobile.setOffline(true); await phone.close();
 const offline=await mobile.newPage(); await offline.goto(base+'/treb/');
 await offline.waitForFunction(()=>!!window.trebViewer);
 await offline.locator('#release').click(); assert.equal(await offline.locator('#range').textContent(),'2177.3 ft');
 // Returning users must migrate out of the previous root-scoped viewer cache.
 serveLegacy=true;
 const migration=await browser.newContext();
 const oldPage=await migration.newPage(); await oldPage.goto(base+'/');
 await oldPage.waitForFunction(async()=>{await navigator.serviceWorker.ready;return !!(await caches.match(location.href));});
 serveLegacy=false;
 await oldPage.evaluate(async()=>{await (await navigator.serviceWorker.getRegistration()).update();});
 await oldPage.waitForFunction(()=>document.title==='First in Fright — Carrboro, NC');
 await oldPage.waitForTimeout(500);
 assert.equal(await oldPage.title(),'First in Fright — Carrboro, NC');
 assert.equal(await oldPage.evaluate(async()=> (await caches.keys()).some(key=>key.startsWith('treb-viewer-'))),false);
 console.log('PASS: homepage desktop/phone; /treb/ viewer scope and offline use; migration from the previous cached homepage without reload loops.');
}finally{await browser.close();await new Promise(done=>server.close(done));}
