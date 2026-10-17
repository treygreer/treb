import {chromium} from 'playwright';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:false});
try {
 const page=await browser.newPage();
 page.on('pageerror',e=>console.log('ERROR',e.message));
 page.on('console',m=>{if(m.type()==='error'||m.type()==='warning')console.log(m.type(),m.text());});
 const start=Date.now();
 await page.goto(pathToFileURL(resolve('site/treb-2026-offline.html')).href);
 try {await page.waitForFunction(()=>!!window.trebViewer,{},{timeout:12000});} catch {}
 console.log(await page.evaluate(()=>({loaded:!!window.trebViewer,details:document.getElementById('details').textContent,status:document.getElementById('status').textContent,webgl:!!document.createElement('canvas').getContext('webgl2')})));
 console.log('Startup',Date.now()-start,'ms');
} finally {await browser.close();}
