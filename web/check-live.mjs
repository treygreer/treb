// Read-only deployment check in an isolated Chrome profile.
import {chromium} from 'playwright';
import {mkdir,writeFile} from 'node:fs/promises';
const url=process.argv[2] && !process.argv[2].startsWith('--') ? process.argv[2] : 'https://treygreer.github.io/treb/';
const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,
  args:['--use-angle=swiftshader','--enable-unsafe-swiftshader']});
try {
 const context=await browser.newContext();
 const hostname=new URL(url).hostname.replace(/^www\./,'');
 const domainAddresses=hostname==='firstinfright.org' ? [`http://${hostname}/`,`http://www.${hostname}/`,`https://${hostname}/`] : [];
 for(const address of domainAddresses) {
  const result=await context.request.get(address,{maxRedirects:0});
  console.log('Domain response:',address,result.status(),'Location:',result.headers().location || '(none)');
 }
 const page=await context.newPage();
 page.on('pageerror',error=>console.log('Script error:',error.message));
 const response=await page.goto(url,{timeout:30000});
 console.log('HTTP:',response.status(),'URL:',page.url());
 console.log('TLS:',await response.securityDetails());
 await page.waitForFunction(()=>!!window.trebViewer,{},{timeout:15000});
 await page.waitForFunction(()=>document.getElementById('status').textContent.startsWith('Offline ready'),{},{timeout:15000});
 console.log(await page.evaluate(()=>({title:document.title,status:document.getElementById('status').textContent,
  viewMode:trebViewer.viewMode,range:trebViewer.data.maxRangeFeet,cache:window.TREB_CACHE,
  slingWeightLb:trebViewer.data.slingWeightLb,armTipHardwareWeightLb:trebViewer.data.armTipHardwareWeightLb,
  sliderMax:document.getElementById('time').max,sliderStep:document.getElementById('time').step})));
 console.log('Page icon links:',await page.locator('link[rel~="icon"]').evaluateAll(links=>links.map(link=>link.href)));
 const favicon=await context.request.get(new URL('favicon.ico',page.url()).href);
 const iconBytes=await favicon.body();
 console.log('favicon.ico HTTP:',favicon.status(),'Type:',favicon.headers()['content-type'],
  'Bytes:',iconBytes.length,'Signature:',iconBytes.subarray(0,12).toString('hex'));
 if(process.argv.includes('--save-icon') && favicon.ok()) {
  await mkdir('src/assets',{recursive:true}); await writeFile('src/assets/favicon.ico',iconBytes);
  console.log('Saved the published icon to src/assets/favicon.ico');
 }
 await context.setOffline(true);
 await page.close();
 const offline=await context.newPage();
 await offline.goto(url,{timeout:15000});
 await offline.waitForFunction(()=>!!window.trebViewer);
 await offline.locator('#release').click();
 console.log('Offline reopen:',await offline.locator('#range').textContent());
} finally {await browser.close();}
