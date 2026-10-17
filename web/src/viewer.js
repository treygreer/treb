import * as THREE from 'three';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';
import {SVGRenderer} from 'three/addons/renderers/SVGRenderer.js';
import data from '../data/viewer.json';

const $ = id => document.getElementById(id);
function unpack(p) {
  const binary=atob(p.base64), bytes=new Uint8Array(binary.length);
  for(let i=0;i<binary.length;i++)bytes[i]=binary.charCodeAt(i);
  return p.dtype === '<u4' ? new Uint32Array(bytes.buffer) : new Float32Array(bytes.buffer);
}
const poses = unpack(data.poses), positions = unpack(data.pumpkin);
const slingPositions = unpack(data.sling), energy = unpack(data.energy);
// Decode forces only when requested, avoiding work and memory during startup.
const forceData = new Map();
let index = 0, selectedTime = data.time[0], playing = false, lastTime = 0, playTime = 0;
const viewport = $('viewport');
const scene = new THREE.Scene();
scene.background = new THREE.Color('#f4f6f9');
const perspectiveCamera = new THREE.PerspectiveCamera(45, 1, .02, 1000);
const flatCamera = new THREE.OrthographicCamera(-10,10,5,-5,.02,1000);
let camera=perspectiveCamera, viewMode='3d';
perspectiveCamera.up.set(0,0,1); flatCamera.up.set(0,0,1);
let renderer, softwareRendering=false;
try {
  renderer = new THREE.WebGLRenderer({antialias:true});
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
} catch (error) {
  softwareRendering=true;
  renderer=new SVGRenderer();
  renderer.setQuality('low');
}
viewport.appendChild(renderer.domElement);
renderer.domElement.setAttribute('aria-label', 'Interactive 3D trebuchet');
renderer.domElement.style.display='block';
renderer.domElement.style.touchAction='none';
let controls;
scene.add(new THREE.HemisphereLight(0xffffff, 0x8c96a3, 2));
const sun = new THREE.DirectionalLight(0xffffff, 1.2);
sun.position.set(-3, -5, 8); scene.add(sun);

const Material=softwareRendering ? THREE.MeshLambertMaterial : THREE.MeshStandardMaterial;
const bodyItems = data.names.map(name => {
  const source = data.meshes[name], verts = unpack(source.vertices), faces = unpack(source.faces);
  const colors = source.colors ? unpack(source.colors) : null;
  const expanded = new Float32Array(faces.length * 3), rgb = new Float32Array(faces.length * 3);
  for (let i = 0; i < faces.length; i++) {
    expanded.set(verts.subarray(faces[i]*3, faces[i]*3+3), i*3);
    const color = colors ? colors.subarray(Math.floor(i/3)*4, Math.floor(i/3)*4+3) : source.color.slice(0,3);
    const linear = new THREE.Color(...color).convertSRGBToLinear();
    rgb.set([linear.r,linear.g,linear.b], i*3);
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(expanded,3));
  geometry.setAttribute('color', new THREE.BufferAttribute(rgb,3));
  geometry.computeVertexNormals();
  const mesh = new THREE.Mesh(geometry, new Material({vertexColors:true, side:THREE.DoubleSide,
    ...(softwareRendering ? {} : {roughness:.85})}));
  mesh.name = name; mesh.matrixAutoUpdate = false;
  scene.add(mesh); return mesh;
});
const pumpkin = new THREE.Mesh(new THREE.SphereGeometry(data.pumpkinRadius, 24, 12),
  new Material({color:0xff6b10}));
scene.add(pumpkin);
function line(color, segments=false) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(6),3));
  const item = new (segments ? THREE.LineSegments : THREE.Line)(geometry, new THREE.LineBasicMaterial({color}));
  item.frustumCulled = false; scene.add(item); return item;
}
const sling = line(0x424955), forces = line(0xdc2635,true), trajectory = line(0xef8834);
const grid = new THREE.GridHelper(18,18,0xb5bbc5,0xd8dde4);
grid.rotation.x = Math.PI/2; grid.position.set(2,0,data.ground); scene.add(grid);
const groundLine=line(0xb5bbc5);
groundLine.geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array([-7,0,data.ground,11,0,data.ground]),3));
groundLine.visible=false;
bodyItems[0].geometry.computeBoundingBox();
const machineBounds=bodyItems[0].geometry.boundingBox;
const viewLow=new THREE.Vector3(Math.min(-3.2,machineBounds.min.x),machineBounds.min.y,data.ground);
const viewHigh=new THREE.Vector3(Math.max(3.2,machineBounds.max.x),machineBounds.max.y,4.6);
const center=viewLow.clone().add(viewHigh).multiplyScalar(.5);
let automaticCamera = true, sideView = false;
let needsRender=true;
const savedViews={};
function createControls(target=center) {
  controls?.dispose();
  controls=new OrbitControls(camera,renderer.domElement);
  controls.target.copy(target);
  controls.enableDamping=true; controls.minDistance=.4; controls.maxDistance=300;
  controls.enableRotate=viewMode==='3d';
  if(viewMode==='2d') {
    controls.mouseButtons.LEFT=THREE.MOUSE.PAN;
    controls.touches.ONE=THREE.TOUCH.PAN;
    controls.touches.TWO=THREE.TOUCH.DOLLY_PAN;
    controls.screenSpacePanning=true;
  }
  controls.addEventListener('change',()=>{needsRender=true;});
  controls.addEventListener('start',()=>{automaticCamera=false;});
  controls.update();
}
createControls();
function flatFrustum() {
  const {width,height}=viewport.getBoundingClientRect(), aspect=width/Math.max(1,height);
  const halfHeight=Math.max((viewHigh.z-viewLow.z)/2,(viewHigh.x-viewLow.x)/2/aspect)*1.08;
  flatCamera.left=-halfHeight*aspect; flatCamera.right=halfHeight*aspect;
  flatCamera.top=halfHeight; flatCamera.bottom=-halfHeight;
  flatCamera.updateProjectionMatrix();
}
function resetCamera(side=false) {
  automaticCamera = true; sideView = side;
  if(viewMode==='2d') {
    flatFrustum(); camera.zoom=1;
    controls.target.copy(center); camera.position.copy(center).add(new THREE.Vector3(0,-40,0));
    camera.updateProjectionMatrix(); createControls(center); needsRender=true;
    return;
  }
  const direction = side ? new THREE.Vector3(0,-1,0) : new THREE.Vector3(.39,-.84,.37).normalize();
  controls.target.copy(center); camera.position.copy(center).add(direction);
  camera.lookAt(center); camera.updateMatrixWorld();
  const tanV=Math.tan(THREE.MathUtils.degToRad(camera.fov/2)), tanH=tanV*camera.aspect;
  let distance=0;
  for(const x of [viewLow.x,viewHigh.x]) for(const y of [viewLow.y,viewHigh.y]) for(const z of [viewLow.z,viewHigh.z]) {
    const point=new THREE.Vector3(x,y,z).applyMatrix4(camera.matrixWorldInverse);
    distance=Math.max(distance,point.z+1+Math.max(Math.abs(point.x)/tanH,Math.abs(point.y)/tanV));
  }
  camera.position.copy(center).addScaledVector(direction,distance*1.08);
  createControls(center);
}
function resize() {
  const {width,height} = viewport.getBoundingClientRect();
  renderer.setSize(width,height);
  perspectiveCamera.aspect=width/height; perspectiveCamera.updateProjectionMatrix();
  flatFrustum();
  if (automaticCamera) resetCamera(sideView);
  needsRender=true;
}
new ResizeObserver(resize).observe(viewport);
function switchDimensions() {
  savedViews[viewMode]={target:controls.target.clone(),automaticCamera,sideView};
  viewMode=viewMode==='3d' ? '2d' : '3d';
  camera=viewMode==='3d' ? perspectiveCamera : flatCamera;
  const saved=savedViews[viewMode];
  createControls(saved?.target);
  if(saved) {automaticCamera=saved.automaticCamera; sideView=saved.sideView;}
  else resetCamera();
  grid.visible=viewMode==='3d'; groundLine.visible=viewMode==='2d';
  $('dimensions').textContent=`2D / 3D: ${viewMode.toUpperCase()}`;
  $('dimensions').setAttribute('aria-pressed',String(viewMode==='2d'));
  $('side').disabled=viewMode==='2d';
  $('camera-help').textContent=viewMode==='2d'
    ? 'Drag to pan · pinch / scroll to zoom'
    : 'Drag to rotate · pinch / scroll to zoom · two fingers / right-drag to pan';
  renderer.domElement.setAttribute('aria-label',`Interactive ${viewMode.toUpperCase()} trebuchet`);
  resize(); needsRender=true;
}

const bars = data.frameNames.map(name => {
  const row = document.createElement('div'); row.className='energy-row';
  const label = document.createElement('span'); label.textContent=name;
  const bar = document.createElement('div'); bar.className='energy-bar';
  const pe=document.createElement('i'), ke=document.createElement('i'); pe.className='pe'; ke.className='ke';
  bar.append(pe,ke); row.append(label,bar); $('energies').append(row);
  return {pe,ke,bar};
});
function setLine(item, values) {
  item.geometry.setAttribute('position',new THREE.BufferAttribute(values,3));
  item.geometry.setDrawRange(0,values.length/3);
}
function sampleAtTime(time) {
  // Search timestamps, not sample numbers: release samples are twenty times
  // denser than the rest of the export. Slider travel must remain linear in time.
  let low=0, high=data.time.length-1;
  while(low<high) {
    const middle=(low+high)>>1;
    if(data.time[middle]<time)low=middle+1;else high=middle;
  }
  return low>0 && time-data.time[low-1]<data.time[low]-time ? low-1 : low;
}
function update(nextIndex, nextTime=data.time[nextIndex]) {
  needsRender=true;
  index = Math.max(0,Math.min(data.time.length-1,nextIndex));
  selectedTime=Math.max(data.time[0],Math.min(data.time.at(-1),nextTime));
  $('time').value=selectedTime;
  $('timestamp').textContent=`${data.time[index].toFixed(4)} s`;
  $('phase').textContent=index<data.releaseIndex ? 'Sling attached' : 'Pumpkin released';
  bodyItems.forEach((mesh,i) => mesh.matrix.set(...poses.subarray((index*bodyItems.length+i)*16,(index*bodyItems.length+i+1)*16)));
  bodyItems[0].visible=$('frame').getAttribute('aria-pressed')==='true';
  pumpkin.position.fromArray(positions,index*3);
  setLine(sling,slingPositions.slice(index*6,index*6+6)); sling.visible=index<data.releaseIndex;
  forces.visible=$('forces').checked;
  if (forces.visible) {
    if(!forceData.has(index))forceData.set(index,unpack(data.forces[index]));
    setLine(forces,forceData.get(index));
  }
  setLine(trajectory,index>=data.releaseIndex ? positions.slice(data.releaseIndex*3,(index+1)*3) : new Float32Array());
  bars.forEach(({pe,ke,bar},i) => {
    const p=energy[(index*bars.length+i)*2], k=energy[(index*bars.length+i)*2+1];
    pe.style.width=`${100*p/data.totalEnergy}%`; ke.style.width=`${100*k/data.totalEnergy}%`;
    bar.title=`Potential: ${p.toFixed(0)} J; kinetic: ${k.toFixed(0)} J`;
  });
  $('range').textContent=`${data.rangeFeet[index].toFixed(1)} ft`;
  $('rangebar').value=Math.max(0,data.rangeFeet[index]);
}
$('time').min=data.time[0];
$('time').max=data.time.at(-1);
$('time').step=data.timeStep;
$('rangebar').max=data.maxRangeFeet;
$('time').addEventListener('input',() => {stop(); const time=Number($('time').value);update(sampleAtTime(time),time);});
$('frame').onclick=() => {$('frame').setAttribute('aria-pressed',String($('frame').getAttribute('aria-pressed')!=='true')); update(index,selectedTime);};
$('forces').onchange=() => update(index,selectedTime);
$('dimensions').onclick=switchDimensions;
$('side').onclick=() => resetCamera(true);
$('reset').onclick=() => resetCamera();
$('release').onclick=() => {stop(); update(data.releaseIndex);};
function stop() {playing=false; $('play').textContent='Play';}
$('play').onclick=() => {
  if (playing) return stop();
  if (index===data.time.length-1) update(0);
  playing=true; playTime=selectedTime; lastTime=performance.now(); $('play').textContent='Pause';
};
function feetAndInches(lengthFeet) {
  // Round to a tenth of an inch before splitting, carrying 12 inches into feet.
  const tenths=Math.round(lengthFeet*120);
  return `${Math.floor(tenths/120)} ft ${(tenths%120)/10} in`;
}
$('details').textContent=`${data.pumpkinWeightLb} lb pumpkin · ${16*data.slingWeightLb} oz sling · ${16*data.armTipHardwareWeightLb} oz added at arm tip · ${feetAndInches(data.slingLengthFeet)} sling · ${data.launchAngle.toFixed(0)}° launch`;
$('release-details').textContent=`Release ${data.releaseTime.toFixed(6)} s · predicted range ${data.maxRangeFeet.toFixed(1)} ft · arm energy ${(data.armEnergyAtRelease/1000).toFixed(2)} kJ`;
function animate(now) {
  if (playing) {
    playTime+=(now-lastTime)/1000 * Number($('speed').value); lastTime=now;
    const next=sampleAtTime(playTime);
    update(next,playTime);
    if(playTime>=data.time.at(-1)) stop();
  }
  controls.update();
  if(needsRender){renderer.render(scene,camera); needsRender=false;}
  requestAnimationFrame(animate);
}
update(0); resize(); requestAnimationFrame(animate);
// Small inspection surface for browser smoke checks; no simulation data mutations.
window.trebViewer={data, get index(){return index;}, get camera(){return camera;}, get viewMode(){return viewMode;}, scene, renderer, update, softwareRendering};
window.trebStartupComplete=true;
$('status').textContent=location.protocol==='file:' ? 'Offline file — everything is included' : 'Viewer loaded';
if(softwareRendering) {
  const note=document.createElement('p');
  note.textContent='Using software 3D rendering because WebGL is unavailable. Camera motion and playback may be slower.';
  $('status').after(note);
}
if ('serviceWorker' in navigator && location.protocol!=='file:' && !window.TREB_DOWNLOAD) {
  navigator.serviceWorker.register('./service-worker.js').then(async () => {
    await navigator.serviceWorker.ready;
    const cache=await caches.open(window.TREB_CACHE);
    if (await cache.match(new URL('./index.html',location.href))) $('status').textContent='Offline ready — you can install this viewer';
  }).catch(error => {$('status').textContent=`Offline installation unavailable: ${error.message}`;});
}
let installPrompt;
window.addEventListener('beforeinstallprompt',event => {event.preventDefault(); installPrompt=event; $('install').hidden=false;});
$('install').onclick=async () => {await installPrompt?.prompt(); $('install').hidden=true;};
window.addEventListener('appinstalled',() => {$('install').hidden=true; $('status').textContent='Installed — ready for offline use';});
