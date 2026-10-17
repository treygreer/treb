import {build} from 'esbuild';
import {readFile,writeFile,mkdir} from 'node:fs/promises';
import {createHash} from 'node:crypto';
const result = await build({entryPoints:['src/viewer.js'],bundle:true,minify:true,format:'iife',target:'es2020',write:false,legalComments:'inline'});
const bundle=result.outputFiles[0].text.replace(/<\/script/gi,'<\\/script');
const template=await readFile('src/template.html','utf8');
const license=await readFile('node_modules/three/LICENSE','utf8');
const favicon=await readFile('src/assets/favicon.ico');
const iconVersion=createHash('sha256').update(favicon).digest('hex').slice(0,12);
const iconUrl=`./favicon.ico?v=${iconVersion}`;
const version=createHash('sha256').update(bundle+template+license).update(favicon).digest('hex').slice(0,12);
const cache=`treb-viewer-${version}`;
function page(download) {
  return template.replace('<!--MANIFEST-->',download?'':'<link rel="manifest" href="./manifest.webmanifest">')
    .replace('<!--FAVICON-->',`<link rel="icon" type="image/x-icon" href="${download ? 'data:image/x-icon;base64,'+favicon.toString('base64') : iconUrl}">`)
    .replace('<!--DOWNLOAD-->',download?'':'<a class="download" href="./treb-2026-offline.html" download>Download Offline HTML</a>')
    .replace('/*SETTINGS*/',`window.TREB_DOWNLOAD=${download};window.TREB_CACHE=${JSON.stringify(cache)};`)
    .replace('/*BUNDLE*/',()=>`/* ${license} */\n${bundle}`);
}
await mkdir('site',{recursive:true});
await writeFile('site/index.html',page(false));
await writeFile('site/treb-2026-offline.html',page(true));
await writeFile('site/favicon.ico',favicon);
await writeFile('site/manifest.webmanifest',JSON.stringify({id:'./',name:'First in Fright 2026',short_name:'Trebuchet',start_url:'./index.html',scope:'./',display:'standalone',background_color:'#eef1f5',theme_color:'#244263',icons:[{src:'icon-192.png',sizes:'192x192',type:'image/png'},{src:'icon-512.png',sizes:'512x512',type:'image/png',purpose:'any maskable'}]},null,2));
await writeFile('site/service-worker.js',`const CACHE=${JSON.stringify(cache)};
const FILES=['./','./index.html','./manifest.webmanifest','./icon-192.png','./icon-512.png','./treb-2026-offline.html','./favicon.ico',${JSON.stringify(iconUrl)}];
self.addEventListener('install',event=>event.waitUntil(caches.open(CACHE).then(cache=>cache.addAll(FILES)).then(()=>self.skipWaiting())));
self.addEventListener('activate',event=>event.waitUntil(caches.keys().then(keys=>Promise.all(keys.filter(key=>key.startsWith('treb-viewer-')&&key!==CACHE).map(key=>caches.delete(key)))).then(()=>self.clients.claim())));
self.addEventListener('fetch',event=>{if(event.request.method!=='GET'||new URL(event.request.url).origin!==self.location.origin)return;event.respondWith(caches.open(CACHE).then(async cache=>(await cache.match(event.request))||fetch(event.request)));});
`);
// esbuild preserves library license comments in the self-contained pages.
await writeFile('site/THREE-LICENSE.txt',license);
console.log(`Built ${version}: two self-contained pages (${(bundle.length/1024/1024).toFixed(2)} MiB each), manifest and offline cache.`);
