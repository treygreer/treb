// Verify the homepage and its separately hosted viewer after publication.
import {chromium} from 'playwright';
import assert from 'node:assert/strict';
const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,
 args:['--use-angle=swiftshader','--enable-unsafe-swiftshader']});
try {
 const context=await browser.newContext();
 const page=await context.newPage();
 const response=await page.goto('https://www.firstinfright.org/');
 assert.equal(response.status(),200);
 assert.equal(await page.title(),'First in Fright — Carrboro, NC');
 assert.equal(await page.locator('.facts strong').first().textContent(),'8 lb');
 const destination=await page.locator('a.button').getAttribute('href');
 assert.equal(destination,'https://treygreer.github.io/treb/');
 await page.locator('a.button').click();
 await page.waitForFunction(()=>!!window.trebViewer);
 await page.waitForFunction(()=>document.getElementById('status').textContent.startsWith('Offline ready'));
 assert.equal(page.url(),destination);
 assert.equal(await page.evaluate(()=>trebViewer.data.slingWeightLb),11/16);
 assert.equal(await page.evaluate(()=>trebViewer.data.pumpkinWeightLb),8);
 assert((await page.locator('#details').textContent()).includes('9 ft 6.5 in sling'));
 assert.equal(await page.evaluate(()=>trebViewer.data.armTipHardwareWeightLb),9/16);
 assert.equal(Number(await page.locator('#time').getAttribute('max')),2);
 for(const time of [.67,.68,.69,.70]) {
  const displayed=await page.evaluate(time=>{
    const slider=document.getElementById('time');slider.value=time;slider.dispatchEvent(new Event('input'));
    return trebViewer.data.time[trebViewer.index];
  },time);
  assert(Math.abs(displayed-time)<=.0011,'Live slider must map uniformly to elapsed seconds');
 }
 await page.locator('#release').click(); assert.equal(await page.locator('#range').textContent(),'2177.3 ft');
 const icon=await context.request.get(new URL(await page.locator('link[rel="icon"]').getAttribute('href'),page.url()).href);
 assert.equal(icon.status(),200);
 await context.setOffline(true); await page.close();
 const offline=await context.newPage(); await offline.goto(destination);
 await offline.waitForFunction(()=>!!window.trebViewer);
 await offline.locator('#release').click(); assert.equal(await offline.locator('#range').textContent(),'2177.3 ft');
 console.log('PASS: live HTTPS homepage, simulation link, updated weights, favicon and offline viewer reopening.');
}finally{await browser.close();}
