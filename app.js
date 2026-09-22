import * as THREE from './three.module.js';

/* ============================================================
   工具函数
============================================================ */
const V = (x,y,z)=>new THREE.Vector3(x,y,z);
const D2R = Math.PI/180;
const S3 = 2*Math.SQRT2/3;          // 0.9428 四面体方向水平分量
const CC = 1.54, CH = 1.09;         // 键长 Å
const VDW = { C:1.70, H:1.20 };
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const lerp=(a,b,t)=>a+(b-a)*t;
const smooth=t=>t*t*(3-2*t);
const angDist=(a,b)=>{let d=Math.abs(a-b)%360;return d>180?360-d:d;};
function csLerp(keys, x){ // 余弦平滑分段插值
  if(x<=keys[0][0]) return keys[0][1];
  for(let i=0;i<keys.length-1;i++){
    const [x0,y0]=keys[i],[x1,y1]=keys[i+1];
    if(x>=x0&&x<=x1){ const u=(x-x0)/(x1-x0); return y0+(y1-y0)*smooth(u); }
  }
  return keys[keys.length-1][1];
}
let tweens=[];
let pumpTimer=null;
function ensurePump(){ // rAF 在后台/嵌入环境可能被暂停，动画期间用定时器兜底
  if(pumpTimer==null)pumpTimer=setInterval(()=>loop(performance.now()),50);
}
function tweenTo(duration, onUpdate, onDone){
  const start=performance.now();
  const tk={start,duration,onUpdate,onDone};
  tweens.push(tk);
  ensurePump();
  return tk;
}
function tickTweens(now){
  tweens=tweens.filter(tk=>{
    let t=(now-tk.start)/tk.duration;
    if(t>=1){ t=1; tk.onUpdate(1); tk.onDone&&tk.onDone(); return false; }
    tk.onUpdate(smooth(t)); return true;
  });
}
function stopTweens(){ tweens=[]; }

/* ============================================================
   分子几何构建 —— 分子 JSON:
   { atoms:[{el, p, tag, label?}], bonds:[[i,j]],
     measures:[{a,b,text,safe}], info:{}, angleData?:[] }
============================================================ */

// ---- 乙烷 ----
function buildEthane(phi){
  const atoms=[], bonds=[], measures=[];
  atoms.push({el:'C',p:V(0,0,0),tag:'C'},{el:'C',p:V(0,0,CC),tag:'C'});
  const fAngs=[90,210,330];
  fAngs.forEach(a=>{ atoms.push({el:'H',p:V(S3*Math.cos(a*D2R)*CH, S3*Math.sin(a*D2R)*CH, -CH/3),tag:'H'}); });
  fAngs.forEach(a0=>{ const a=a0+phi;
    atoms.push({el:'H',p:V(S3*Math.cos(a*D2R)*CH, S3*Math.sin(a*D2R)*CH, CC+CH/3),tag:'H'}); });
  bonds.push([0,1]);
  for(let i=0;i<3;i++){ bonds.push([0,2+i],[1,5+i]); }
  // 最近的前-后 H 距离
  let hh=99, ai=-1, bi=-1;
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){
    const d=atoms[2+i].p.distanceTo(atoms[5+j].p);
    if(d<hh){hh=d;ai=2+i;bi=5+j;}
  }
  measures.push({a:ai,b:bi,text:(hh*100).toFixed(0)+' pm',safe:hh>=2.40});
  return {atoms,bonds,measures,info:{hh}, front:0, back:1};
}

// ---- 甲基（在中心碳 M 上，键轴 u 指向所连中心碳）----
function attachMethyl(atoms, bonds, M, u, chainTag){
  const e1 = Math.abs(u.y)>0.9 ? V(1,0,0) : V(0,1,0).sub(u.clone().multiplyScalar(u.y)).normalize();
  const e2 = new THREE.Vector3().crossVectors(u,e1).normalize();
  const mc = atoms.length;
  atoms.push({el:'C',p:M,tag:chainTag==='tb'?'tb':'me',label:chainTag==='tb'?'CH₃':'Me'});
  for(let i=0;i<3;i++){
    const a=i*120*D2R;
    const dir=u.clone().multiplyScalar(-1/3)
      .add(e1.clone().multiplyScalar(S3*Math.cos(a)))
      .add(e2.clone().multiplyScalar(S3*Math.sin(a))).normalize();
    atoms.push({el:'H',p:M.clone().add(dir.multiplyScalar(CH)),tag:chainTag==='tb'?'tbH':'meH'});
    bonds.push([mc,atoms.length-1]);
  }
  return mc;
}
// ---- 叔丁基 ----
function attachTBu(atoms,bonds,S,u){
  const e1 = Math.abs(u.y)>0.9 ? V(1,0,0) : V(0,1,0).sub(u.clone().multiplyScalar(u.y)).normalize();
  const e2 = new THREE.Vector3().crossVectors(u,e1).normalize();
  const cen=atoms.length;
  atoms.push({el:'C',p:S,tag:'tb',label:'t-Bu'});
  for(let i=0;i<3;i++){
    const a=(i*120+30)*D2R;
    const dir=u.clone().multiplyScalar(-1/3)
      .add(e1.clone().multiplyScalar(S3*Math.cos(a)))
      .add(e2.clone().multiplyScalar(S3*Math.sin(a))).normalize();
    const Mp=S.clone().add(dir.clone().multiplyScalar(CC));
    const mc=attachMethyl(atoms,bonds,Mp,dir,'tb');
    bonds.push([cen,mc]);
  }
  return cen;
}

// ---- 丁烷（theta = 两甲基二面角）----
function buildButane(theta){
  const atoms=[],bonds=[],measures=[];
  atoms.push({el:'C',p:V(0,0,0),tag:'C',label:'C2'},{el:'C',p:V(0,0,CC),tag:'C',label:'C3'});
  const f=[{a:90,t:'me'},{a:210,t:'H'},{a:330,t:'H'}];
  const b=[{a:90+theta,t:'me'},{a:210+theta,t:'H'},{a:330+theta,t:'H'}];
  let meF=-1,meB=-1;
  f.forEach(s=>{
    const d=V(S3*Math.cos(s.a*D2R),S3*Math.sin(s.a*D2R),-1/3);
    if(s.t==='me'){
      const M=d.clone().multiplyScalar(CC);
      meF=attachMethyl(atoms,bonds,M,d.clone().negate(),'me');
      bonds.push([0,meF]);
    } else { atoms.push({el:'H',p:d.multiplyScalar(CH),tag:'H'}); bonds.push([0,atoms.length-1]); }
  });
  b.forEach(s=>{
    const d=V(S3*Math.cos(s.a*D2R),S3*Math.sin(s.a*D2R), 1/3);
    if(s.t==='me'){
      const M=V(0,0,CC).add(d.clone().multiplyScalar(CC));
      meB=attachMethyl(atoms,bonds,M,d.clone().negate(),'me');
      bonds.push([1,meB]);
    } else { atoms.push({el:'H',p:V(0,0,CC).add(d.multiplyScalar(CH)),tag:'H'}); bonds.push([1,atoms.length-1]); }
  });
  bonds.push([0,1]);
  const mm=atoms[meF].p.distanceTo(atoms[meB].p);
  measures.push({a:meF,b:meB,text:(mm*100).toFixed(0)+' pm',safe:mm>=3.40});
  return {atoms,bonds,measures,info:{mm},front:0,back:1};
}

// ---- 环上碳的两根 H（局部四面体）----
function ringHDirs(p,n1,n2){
  const u1=n1.clone().sub(p).normalize(), u2=n2.clone().sub(p).normalize();
  const bis=u1.clone().add(u2);
  if(bis.lengthSq()<1e-6) return [V(0,0,1),V(0,0,-1)];
  bis.normalize();
  let n=new THREE.Vector3().crossVectors(u1,u2);
  if(n.lengthSq()<1e-6) n=V(0,0,1); else n.normalize();
  const c1=1/Math.sqrt(3), c2=Math.sqrt(2/3);
  const h1=bis.clone().multiplyScalar(-c1).add(n.clone().multiplyScalar(c2));
  const h2=bis.clone().multiplyScalar(-c1).add(n.clone().multiplyScalar(-c2));
  return [h1.normalize(),h2.normalize()];
}

// ---- 环烷烃（3/4/5/6椅）----
function buildRing(n){
  const atoms=[],bonds=[],cs=[];
  let angleData=[], avgAngle=0, minAngle=180;
  if(n===3){
    const R=CC/Math.sqrt(3);
    [90,210,330].forEach(a=>cs.push(V(R*Math.cos(a*D2R),R*Math.sin(a*D2R),0)));
  }else if(n===4){
    const R=CC/Math.SQRT2;
    [45,135,225,315].forEach((a,k)=>cs.push(V(R*Math.cos(a*D2R),R*Math.sin(a*D2R),k%2? .24:-.24)));
  }else if(n===5){
    const R=1.272;
    [90,162,234,306,18].forEach((a,k)=>cs.push(V(R*Math.cos(a*D2R),R*Math.sin(a*D2R),k===0?.55:0)));
  }else{
    chairFrame(0).forEach(p=>cs.push(p));
  }
  cs.forEach(p=>atoms.push({el:'C',p:p.clone(),tag:'C'}));
  for(let k=0;k<n;k++){
    bonds.push([k,(k+1)%n]);
    const p=cs[k], p1=cs[(k+1)%n], p2=cs[(k-1+n)%n];
    let dirs;
    if(n===3){ // 平面内向外张开
      const out=p.clone().negate().normalize();
      const base=Math.atan2(out.y,out.x);
      dirs=[57,-57].map(d=>V(Math.cos(base+d*D2R),Math.sin(base+d*D2R),0));
    }else dirs=ringHDirs(p,p1,p2);
    dirs.forEach(d=>atoms.push({el:'H',p:p.clone().add(d.multiplyScalar(CH)),tag:'H'}));
    const ang=p1.clone().sub(p).angleTo(p2.clone().sub(p))/D2R;
    avgAngle+=ang/n; minAngle=Math.min(minAngle,ang);
    angleData.push({at:k,n1:(k+1)%n,n2:(k-1+n)%n,v:ang});
  }
  for(let k=0;k<n;k++){ bonds.push([k,n+2*k],[k,n+2*k+1]); }
  return {atoms,bonds,measures:[],info:{avgAngle,minAngle},angleData};
}

/* ---- 环己烷椅式/船式构象关键帧 ----
   顶点 k=0..5，基础方位角 k*60°，r 为环半径，z 为上下位移，tw 为扭转角 */
const CHAIR_KFS=[
  {t:0.00, r:1.452, z:[ .257,-.257, .257,-.257, .257,-.257], tw:[0,0,0,0,0,0],         name:'椅式'},
  {t:0.17, r:1.430, z:[ .55 , .10 ,-.08 ,-.12 ,-.08 , .10 ], tw:[0,1.5,-1.5,0,-1.5,1.5],name:'半椅式'},
  {t:0.34, r:1.420, z:[ .58 , .25 ,-.25 , .58 ,-.25 , .25 ], tw:[7,3,-3,-7,-3,3],       name:'扭船式'},
  {t:0.50, r:1.406, z:[ .63 , 0   , 0   , .63 , 0   , 0   ], tw:[0,0,0,0,0,0],         name:'船式'},
  {t:0.66, r:1.420, z:[ .58 , .25 ,-.25 , .58 ,-.25 , .25 ], tw:[-7,-3,3,7,3,-3],      name:'扭船式'},
  {t:0.83, r:1.430, z:[-.10 ,-.08 , .10 ,-.55 , .10 ,-.08 ], tw:[0,-1.5,1.5,0,1.5,-1.5],name:'半椅式'},
  {t:1.00, r:1.452, z:[-.257, .257,-.257, .257,-.257, .257], tw:[0,0,0,0,0,0],         name:'椅式′'},
];
const CHAIR_E_KEYS=[[0,0],[17,50.6],[34,23],[50,30],[66,23],[83,50.6],[100,0]];
const CHAIR_NAME_KEYS=CHAIR_KFS.map(k=>[k.t*100,k.name]);

function chairFrame(t){ // 返回 6 个碳坐标
  let kf0=CHAIR_KFS[0],kf1=CHAIR_KFS[CHAIR_KFS.length-1];
  for(let i=0;i<CHAIR_KFS.length-1;i++){
    if(t>=CHAIR_KFS[i].t&&t<=CHAIR_KFS[i+1].t){kf0=CHAIR_KFS[i];kf1=CHAIR_KFS[i+1];break;}
  }
  const u=smooth(clamp((t-kf0.t)/(kf1.t-kf0.t),0,1));
  const r=lerp(kf0.r,kf1.r,u);
  return Array.from({length:6},(_,k)=>{
    const a=(k*60+lerp(kf0.tw[k],kf1.tw[k],u))*D2R;
    return V(r*Math.cos(a),r*Math.sin(a),lerp(kf0.z[k],kf1.z[k],u));
  });
}

function buildChairMorph(t){
  const cs=chairFrame(t);
  const atoms=[],bonds=[];
  cs.forEach(p=>atoms.push({el:'C',p:p.clone(),tag:'C'}));
  for(let k=0;k<6;k++) bonds.push([k,(k+1)%6]);
  // H：在船式区段，C0/C3 给出旗杆氢
  const boatness=clamp(1-Math.abs(t-0.5)/0.28,0,1);
  for(let k=0;k<6;k++){
    const p=cs[k], p1=cs[(k+1)%6], p2=cs[(k+5)%6];
    let dirs=ringHDirs(p,p1,p2);
    if(boatness>0.05&&(k===0||k===3)){
      const a0=k*60*D2R;
      const radial=V(Math.cos(a0),Math.sin(a0),0);
      const flag=radial.clone().multiplyScalar(-.45).add(V(0,0,.893)).normalize();
      const other=flag.clone().negate().sub(p1.clone().sub(p).normalize()).sub(p2.clone().sub(p).normalize());
      if(other.lengthSq()<1e-6) other.set(0,0,-1); else other.normalize();
      // 随 boatness 混合，保证椅式端仍是四面体方向
      const fIdx=dirs.reduce((best,d,i)=>d.z>dirs[best].z?i:best,0);
      const oIdx=1-fIdx;
      const df=dirs[fIdx].clone().lerp(flag,boatness);
      const dO=dirs[oIdx].clone().lerp(other,boatness);
      if(df.lengthSq()<1e-6)df.copy(flag); if(dO.lengthSq()<1e-6)dO.copy(other);
      dirs=[df.normalize(),dO.normalize()];
    }
    dirs.forEach(d=>atoms.push({el:'H',p:p.clone().add(d.clone().multiplyScalar(CH)),tag:'H',_carbon:k}));
  }
  for(let k=0;k<6;k++){ bonds.push([k,6+2*k],[k,6+2*k+1]); }
  // 旗杆 H 距离：C0 与 C3 上最靠近的一对 H
  let best=99,bi=-1,bj=-1;
  for(let i=0;i<2;i++)for(let j=0;j<2;j++){
    const a=6+0*2+i, b=6+3*2+j, d=atoms[a].p.distanceTo(atoms[b].p);
    if(d<best){best=d;bi=a;bj=b;}
  }
  if(boatness>0.35&&bi>=0){ atoms[bi].tag='flag'; atoms[bj].tag='flag'; }
  const measures=boatness>0.12&&bi>=0?[{a:bi,b:bj,text:(best*100).toFixed(0)+' pm',safe:best>=2.40}]:[];
  return {atoms,bonds,measures,info:{flag:best,boatness},front:1,back:2};
}

/* ---- 分析型椅式（a/e 键严格准确），flip=0/1 ----
   偶数碳在上层（a 朝上），奇数碳下层（a 朝下）；flip 后互换 */
function analyticChair(flip){
  const r=1.452,h=.257;
  const C=[];
  for(let k=0;k<6;k++){
    const a=k*60*D2R;
    let z=(k%2? -h:h); if(flip) z=-z;
    C.push(V(r*Math.cos(a),r*Math.sin(a),z));
  }
  // 每碳两个键方向: {a, e, up}
  const dirs=[];
  for(let k=0;k<6;k++){
    const a=k*60*D2R, axUp=k%2===0;
    const upAx=axUp?!flip:flip;   // 该碳直立键是否朝上
    const aDir=V(0,0,upAx?1:-1);
    const eDir=V(S3*Math.cos(a),S3*Math.sin(a),upAx?-1/3:1/3);
    dirs.push({a:aDir,e:eDir,aUp:upAx,eUp:!upAx});
  }
  return {C,dirs};
}

/* 取代椅式
   subs:[{k, kind:'Me'|'tBu'|'H', face:'up'|'down'}], t 为翻环进度 0..1 */
function buildSubChair(subs, t){
  const f0=analyticChair(false), f1=analyticChair(true);
  const cs=[], dirsEnd0=[], dirsEnd1=[];
  for(let k=0;k<6;k++){
    cs.push(f0.C[k].clone().lerp(f1.C[k],smooth(t)));
    dirsEnd0.push(f0.dirs[k]); dirsEnd1.push(f1.dirs[k]);
  }
  const atoms=[],bonds=[];
  cs.forEach(p=>atoms.push({el:'C',p:p.clone(),tag:'C'}));
  for(let k=0;k<6;k++) bonds.push([k,(k+1)%6]);
  const subMap={}; subs.forEach(s=>subMap[s.k]=s);
  const contacts=[];
  for(let k=0;k<6;k++){
    const p=cs[k], p1=cs[(k+1)%6], p2=cs[(k+5)%6];
    const s=subMap[k];
    // 该碳“非环键”在 t=0 和 t=1 的方向
    const slot=(face)=>{
      const d0init = face==='up' ? (dirsEnd0[k].aUp?dirsEnd0[k].a:dirsEnd0[k].e)
                                  : (dirsEnd0[k].aUp?dirsEnd0[k].e:dirsEnd0[k].a);
      const d1end  = face==='up' ? (dirsEnd1[k].aUp?dirsEnd1[k].a:dirsEnd1[k].e)
                                  : (dirsEnd1[k].aUp?dirsEnd1[k].e:dirsEnd1[k].a);
      // 球面插值
      const d=d0init.clone().lerp(d1end,smooth(t));
      if(d.lengthSq()<1e-6) d.copy(d0init);
      return {d:d.normalize(), startType: face==='up' ? (dirsEnd0[k].aUp?'a':'e') : (dirsEnd0[k].aUp?'e':'a'),
                        endType: face==='up' ? (dirsEnd1[k].aUp?'a':'e') : (dirsEnd1[k].aUp?'e':'a')};
    };
    if(s&&s.kind!=='H'){
      const {d,startType,endType}=slot(s.face);
      const typeNow = t<.5?startType:endType;
      const u=d.clone().negate();
      if(s.kind==='Me'){
        const M=p.clone().add(d.clone().multiplyScalar(CC));
        const mc=attachMethyl(atoms,bonds,M,u,'me');
        bonds.push([k,mc]);
        atoms[mc].slot=typeNow; atoms[mc]._sub=true;
      }else{
        const S=p.clone().add(d.clone().multiplyScalar(CC));
        const cen=attachTBu(atoms,bonds,S,u);
        bonds.push([k,cen]);
        atoms[cen].slot=typeNow; atoms[cen]._sub=true;
      }
      // 补一个环 H（取与取代方向偏差最大的四面体方向）
      const dd=ringHDirs(p,p1,p2).sort((x,y)=>x.angleTo(d)-y.angleTo(d)).pop();
      atoms.push({el:'H',p:p.clone().add(dd.multiplyScalar(CH)),tag:'H',_carbon:k});
      bonds.push([k,atoms.length-1]);
    }else{
      let [dA,dE]=[
        slot('up').d, slot('down').d
      ];
      const up0=dirsEnd0[k].aUp;
      atoms.push({el:'H',p:p.clone().add(dA.clone().multiplyScalar(CH)),tag:(t<.5?(up0?'ax':'eq'):(up0?'eq':'ax')),_carbon:k,face:'up'});
      atoms.push({el:'H',p:p.clone().add(dE.clone().multiplyScalar(CH)),tag:(t<.5?(up0?'eq':'ax'):(up0?'ax':'eq')),_carbon:k,face:'down'});
      bonds.push([k,atoms.length-2],[k,atoms.length-1]);
    }
  }
  return {atoms,bonds,measures:[],info:{}};
}

/* ============================================================
   Three.js 分子查看器
============================================================ */
const COLORS={
  C:0xaeb9c9, H:0xeef3fb, me:0x5aa7ff, meH:0xb9dcff, tb:0x3fd4c0, tbH:0xbef5ee,
  ax:0x6ea8ff, eq:0x46d6c8, flag:0xff8fa0
};
const sphereGeo={C:new THREE.SphereGeometry(.31,28,22), H:new THREE.SphereGeometry(.19,20,16)};
const vdwGeo={C:new THREE.SphereGeometry(VDW.C,26,20), H:new THREE.SphereGeometry(VDW.H,22,16)};
const bondGeo=new THREE.CylinderGeometry(1,1,1,14);

class Viewer{
  constructor(wrapId, cam0){
    this.wrap=document.getElementById(wrapId);
    this.canvas=this.wrap.querySelector('.mol-canvas');
    this.ov=this.wrap.querySelector('.ov-canvas');
    this.ctx=this.ov.getContext('2d');
    this.renderer=new THREE.WebGLRenderer({canvas:this.canvas,antialias:true,alpha:true});
    this.renderer.setPixelRatio(Math.min(devicePixelRatio,2));
    this.scene=new THREE.Scene();
    this.camera=new THREE.PerspectiveCamera(42,1,0.1,100);
    this.cam0=cam0||{r:7.2,phi:1.18,theta:0.85,tx:0,ty:0,tz:0.1};
    this.cam={...this.cam0};
    this.opts={vdw:true,vdwScale:.6,rotate:true,labels:false,angleLabel:false,flags:false,aeColor:true};
    this._lastInteract=0;
    this.group=new THREE.Group(); this.scene.add(this.group);
    this.lg=new THREE.Group(); this.scene.add(this.lg);
    const amb=new THREE.HemisphereLight(0xbfd8ff,0x141c2c,1.05); this.scene.add(amb);
    const dl=new THREE.DirectionalLight(0xffffff,1.7); dl.position.set(4,7,5); this.scene.add(dl);
    const dl2=new THREE.DirectionalLight(0x6f9fff,.55); dl2.position.set(-6,-3,-4); this.scene.add(dl2);
    this.mol=null; this.clashes=[];
    this._bind();
    new ResizeObserver(()=>this._resize()).observe(this.wrap);
    this._resize();
  }
  resetView(){ stopTweens(); const a=this.cam0,b=this.cam;
    tweenTo(420,t=>{b.r=lerp(b.r,a.r,t);b.phi=lerp(b.phi,a.phi,t);b.theta=lerp(b.theta,a.theta,t);b.tz=lerp(b.tz,a.tz,t);}); }
  _resize(){
    const w=this.wrap.clientWidth,h=this.wrap.clientHeight;
    if(!w)return;
    this.renderer.setSize(w,h,false);
    this.ov.width=w*devicePixelRatio; this.ov.height=h*devicePixelRatio;
    this.ov.style.width=w+'px'; this.ov.style.height=h+'px';
    this.ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
    this.camera.aspect=w/h; this.camera.updateProjectionMatrix();
    this.w=w; this.h=h;
  }
  _bind(){
    let drag=false,lx=0,ly=0,pinch=0;
    const down=e=>{drag=true;lx=e.clientX;ly=e.clientY;this._lastInteract=performance.now();};
    this.canvas.addEventListener('pointerdown',e=>{this.canvas.setPointerCapture(e.pointerId);down(e);});
    this.canvas.addEventListener('pointermove',e=>{
      if(!drag)return;
      const dx=e.clientX-lx,dy=e.clientY-ly;lx=e.clientX;ly=e.clientY;
      if(e.pointerType==='touch'&&this._d2)return;
      this.cam.theta-=dx*.006; this.cam.phi=clamp(this.cam.phi-dy*.006,.15,Math.PI-.15);
      this._lastInteract=performance.now();
    });
    const up=()=>{if(drag){this._lastInteract=performance.now();}drag=false;};
    window.addEventListener('pointerup',up);
    this.canvas.addEventListener('wheel',e=>{e.preventDefault();
      this.cam.r=clamp(this.cam.r*(1+e.deltaY*.0011),3.2,20);this._lastInteract=performance.now();},{passive:false});
    this.canvas.addEventListener('touchstart',e=>{
      if(e.touches.length===2){const dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;pinch=Math.hypot(dx,dy);this._d2=true;drag=false;this._lastInteract=performance.now();}
    },{passive:true});
    this.canvas.addEventListener('touchmove',e=>{
      if(e.touches.length===2){const dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;
        const d=Math.hypot(dx,dy);this.cam.r=clamp(this.cam.r*pinch/d,3.2,20);pinch=d;this._lastInteract=performance.now();}
    },{passive:true});
    this.canvas.addEventListener('touchend',()=>this._d2=false);
  }
  clearGroup(g){ while(g.children.length){ const c=g.children.pop(); g.remove(c);
    if(c.geometry&&!Object.values(sphereGeo).includes(c.geometry)&&!Object.values(vdwGeo).includes(c.geometry)&&c.geometry!==bondGeo) c.geometry.dispose();
    if(c.material){ if(Array.isArray(c.material))c.material.forEach(m=>m.dispose()); else c.material.dispose(); }
  }}
  setMolecule(mol){
    this.mol=mol;
    this.clearGroup(this.group); this.clearGroup(this.lg);
    // 非键近距扫描：
    //  键连路径 ≤2 键（1-2/1-3）一律排除；
    //  路径恰好 3 键（1-4，如环骨架 C…C、重叠式 H-C-C-H）只计 H…H（扭转重叠）；
    //  路径 ≥4 键（旗杆 H、1,3-二直立、甲基互斥等真正的空间位阻）全部计入。
    const adj=mol.atoms.map(()=>[]);
    mol.bonds.forEach(([a,b])=>{adj[a].push(b);adj[b].push(a);});
    const bondDist=(s,target)=>{ // BFS，最多 4 层
      let frontier=[s],seen=new Set([s]);
      for(let depth=1;depth<=4;depth++){
        const nxt=[];
        for(const x of frontier)for(const y of adj[x]){
          if(y===target)return depth;
          if(!seen.has(y)){seen.add(y);nxt.push(y);}
        }
        frontier=nxt;
      }
      return 5;
    };
    const atomSev=new Map(); this.clashes=[];
    for(let i=0;i<mol.atoms.length;i++)for(let j=i+1;j<mol.atoms.length;j++){
      const bd=bondDist(i,j);
      if(bd<=2)continue;
      if(bd===3&&!(mol.atoms[i].el==='H'&&mol.atoms[j].el==='H'))continue;
      const d=mol.atoms[i].p.distanceTo(mol.atoms[j].p);
      const lim=VDW[mol.atoms[i].el]+VDW[mol.atoms[j].el];
      if(d<lim-.01){
        const sev=clamp((lim-d)/.9,0,1);
        this.clashes.push({i,j,d,sev});
        atomSev.set(i,Math.max(atomSev.get(i)||0,sev));
        atomSev.set(j,Math.max(atomSev.get(j)||0,sev));
      }
    }
    this.atomSev=atomSev;
    // 键
    mol.bonds.forEach(([a,b])=>{
      const pa=mol.atoms[a].p,pb=mol.atoms[b].p;
      const heavy=mol.atoms[a].el==='C'&&mol.atoms[b].el==='C';
      const m=new THREE.Mesh(bondGeo,new THREE.MeshPhongMaterial({color:0x8fa3bd,shininess:40}));
      const mid=pa.clone().add(pb).multiplyScalar(.5), dir=pb.clone().sub(pa);
      m.position.copy(mid); m.quaternion.setFromUnitVectors(V(0,1,0),dir.clone().normalize());
      m.scale.set(heavy?.075:.052,dir.length(),heavy?.075:.052);
      this.group.add(m);
    });
    // 原子
    mol.atoms.forEach((at,i)=>{
      let tag=at.tag;
      if(!this.opts.aeColor&&(tag==='ax'||tag==='eq')) tag='H';
      const col=COLORS[tag]??COLORS[at.el];
      const m=new THREE.Mesh(sphereGeo[at.el],new THREE.MeshPhongMaterial({color:col,shininess:55,
        specular:0x33445c, emissive:tag==='flag'&&this.opts.flags?0x55101b:0x000000}));
      m.position.copy(at.p); this.group.add(m);
      if(this.opts.vdw && this.opts.vdwScale>0){
        const sev=atomSev.get(i)||0;
        const strong=tag==='flag'||sev>=.45;          // 强空间冲突才染红
        const isBig=(tag==='me'||tag==='tb');
        const opScale=Math.min(1,this.opts.vdwScale+.18);   // 缩小时保留可见度
        const vm=new THREE.Mesh(vdwGeo[at.el],new THREE.MeshPhongMaterial({
          color:strong?0xff4d66:col, transparent:true,
          opacity:(strong?.22:(isBig?.16:.10))*opScale,
          depthWrite:false, shininess:90}));
        vm.position.copy(at.p);
        if(this.opts.vdwScale!==1) vm.scale.setScalar(this.opts.vdwScale);
        this.group.add(vm);
      }
      if(this.opts.labels&&at.label){
        const sp=makeSprite(at.label, tag==='me'?'#8fc2ff':tag==='tb'?'#7feadd':'#cfe0f5');
        sp.position.copy(at.p).add(V(0,.42,0)); sp.scale.set(.62,.31,1); this.lg.add(sp);
      }
      if(this.opts.aeColor&&(tag==='ax'||tag==='eq')){
        const sp=makeSprite(tag==='ax'?'a':'e', tag==='ax'?'#8fbcff':'#76e7d6');
        sp.position.copy(at.p).add(V(0,.30,0)); sp.scale.set(.34,.27,1); this.lg.add(sp);
      }
      if(this.opts.flags&&tag==='flag'){
        const sp=makeSprite('旗杆H','#ffb3bd');
        sp.position.copy(at.p).add(V(0,.34,0)); sp.scale.set(.62,.30,1); this.lg.add(sp);
      }
    });
    this.render();   // 立即上屏一帧（后台 rAF 暂停时也保证交互即时可见）
  }
  project(p){
    const v=p.clone().project(this.camera);
    return {x:(v.x*.5+.5)*this.w, y:(-v.y*.5+.5)*this.h, behind:v.z>1};
  }
  render(){
    if(!this.wrap.offsetParent) return;
    const now=performance.now();
    if(this.opts.rotate && now-this._lastInteract>1500) this.cam.theta+=.0032;
    const {r,phi,theta,tx,ty,tz}=this.cam;
    this.camera.position.set(tx+r*Math.sin(phi)*Math.sin(theta), ty+r*Math.cos(phi), tz+r*Math.sin(phi)*Math.cos(theta));
    this.camera.lookAt(tx,ty,tz);
    this.renderer.render(this.scene,this.camera);
    // overlay
    const ctx=this.ctx; ctx.clearRect(0,0,this.w,this.h);
    if(!this.mol)return;
    const drawDash=(a,b,color,width,gap)=>{
      ctx.save();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(gap);
      ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke();ctx.restore();
    };
    // 冲突接触：强冲突（旗杆氢等）红色醒目，轻微重叠接触淡橙色提示扭转张力
    this.clashes.forEach(c=>{
      if(c.sev<.05)return;
      const strong=c.sev>=.45;
      const A=this.project(this.mol.atoms[c.i].p),B=this.project(this.mol.atoms[c.j].p);
      if(A.behind||B.behind)return;
      if(strong)drawDash(A,B,`rgba(255,80,105,${.5+.5*c.sev})`,1.6+1.6*c.sev,[5,4]);
      else drawDash(A,B,`rgba(255,165,95,${.14+.28*c.sev})`,1,[4,4]);
    });
    // 测量距离
    this.mol.measures.forEach(m=>{
      const pa=this.mol.atoms[m.a]?.p,pb=this.mol.atoms[m.b]?.p;
      if(!pa||!pb)return;
      const A=this.project(pa),B=this.project(pb);
      if(A.behind||B.behind)return;
      const col=m.safe?'rgba(74,222,174,.9)':'rgba(255,99,120,.95)';
      drawDash(A,B,col,2,[7,4]);
      const mx=(A.x+B.x)/2,my=(A.y+B.y)/2;
      labelPill(ctx,mx,my-9,m.text,col);
    });
    // 键角标注
    if(this.opts.angleLabel&&this.mol.angleData){
      this.mol.angleData.forEach(d=>{
        const C=this.project(this.mol.atoms[d.at].p);
        const N1=this.project(this.mol.atoms[d.n1].p);
        const N2=this.project(this.mol.atoms[d.n2].p);
        if(C.behind)return;
        let a1=Math.atan2(N1.y-C.y,N1.x-C.x), a2=Math.atan2(N2.y-C.y,N2.x-C.x);
        const hot=d.v<100;
        ctx.strokeStyle=hot?'rgba(255,99,120,.8)':'rgba(140,180,255,.55)';
        ctx.lineWidth=1.4;ctx.beginPath();ctx.arc(C.x,C.y,20,a1,a2);ctx.stroke();
        const aa=(a1+a2)/2;
        labelPill(ctx,C.x+Math.cos(aa)*34,C.y+Math.sin(aa)*34,Math.round(d.v)+'°',
          hot?'rgba(255,99,120,.95)':'rgba(150,190,255,.9)');
      });
    }
  }
}

function labelPill(ctx,x,y,text,color){
  ctx.save();ctx.font='600 11.5px -apple-system,sans-serif';
  const w=ctx.measureText(text).width+14;
  ctx.fillStyle='rgba(13,20,32,.82)';ctx.strokeStyle=color;ctx.lineWidth=1;
  roundRect(ctx,x-w/2,y-10,w,20,6);ctx.fill();ctx.stroke();
  ctx.fillStyle=color;ctx.textAlign='center';ctx.textBaseline='middle';
  ctx.fillText(text,x,y+1);ctx.restore();
}
function roundRect(ctx,x,y,w,h,r){ctx.beginPath();ctx.moveTo(x+r,y);ctx.arcTo(x+w,y,x+w,y+h,r);
  ctx.arcTo(x+w,y+h,x,y+h,r);ctx.arcTo(x,y+h,x,y,r);ctx.arcTo(x,y,x+w,y,r);ctx.closePath();}

const spriteCache=new Map();
function makeSprite(text,color){
  const key=text+color;
  if(spriteCache.has(key))return spriteCache.get(key).clone();
  const c=document.createElement('canvas');c.width=256;c.height=128;
  const x=c.getContext('2d');
  x.font='600 52px -apple-system,"PingFang SC",sans-serif';x.textAlign='center';x.textBaseline='middle';
  x.shadowColor='rgba(0,0,0,.8)';x.shadowBlur=8;
  x.fillStyle=color;x.fillText(text,128,66);
  const tex=new THREE.CanvasTexture(c);tex.anisotropy=4;
  const sp=new THREE.Sprite(new THREE.SpriteMaterial({map:tex,transparent:true,depthTest:false}));
  spriteCache.set(key,sp);return sp.clone();
}

/* ============================================================
   Newman 投影式
============================================================ */
function drawNewman(canvas, mol, frontId, backId){
  const dpr=Math.min(devicePixelRatio,2);
  const W=canvas.clientWidth, H=canvas.clientHeight;
  canvas.width=W*dpr;canvas.height=H*dpr;
  const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.clearRect(0,0,W,H);
  const cx=W/2,cy=H/2,L=Math.min(W,H)*.34,gap=12;
  const F=mol.atoms[frontId].p,B=mol.atoms[backId].p;
  const axis=B.clone().sub(F).normalize();
  let e1=V(0,1,0).sub(axis.clone().multiplyScalar(axis.y));
  if(e1.lengthSq()<1e-4)e1=V(1,0,0);
  e1.normalize();
  const e2=new THREE.Vector3().crossVectors(axis,e1).normalize();
  const bonded=new Map();
  mol.bonds.forEach(([a,b])=>{
    if(!bonded.has(a))bonded.set(a,[]); if(!bonded.has(b))bonded.set(b,[]);
    bonded.get(a).push(b);bonded.get(b).push(a);
  });
  const subOf=id=>bonded.get(id).filter(x=>x!==frontId&&x!==backId).map(sid=>{
    const dir=mol.atoms[sid].p.clone().sub(id===frontId?F:B);
    return {x:dir.dot(e1),y:dir.dot(e2),at:mol.atoms[sid]};
  });
  const tagCol=t=>t==='me'||t==='tb'?'#5aa7ff':t==='meH'||t==='tbH'?'#9cc7ff':'#d9e4f4';
  const tagName=at=>{
    if(at.tag==='me')return 'CH₃';
    if(at.tag==='tb')return 'tBu';
    if(at.el==='C')return 'CH₂';
    return 'H';
  };
  // 后碳：圆 + 三条线
  ctx.strokeStyle='#8aa1c0';ctx.lineWidth=2;
  ctx.beginPath();ctx.arc(cx,cy,gap,0,Math.PI*2);ctx.stroke();
  subOf(backId).forEach(s=>{
    const a=Math.atan2(s.y,s.x);
    const x0=cx+Math.cos(a)*(gap+2),y0=cy-Math.sin(a)*(gap+2);
    const x1=cx+Math.cos(a)*L,y1=cy-Math.sin(a)*L;
    ctx.strokeStyle=tagCol(s.at.tag);ctx.lineWidth=2;
    ctx.beginPath();ctx.moveTo(x0,y0);ctx.lineTo(x1,y1);ctx.stroke();
    ctx.fillStyle=tagCol(s.at.tag);ctx.font='600 13px -apple-system,sans-serif';
    ctx.textAlign='center';ctx.textBaseline='middle';
    ctx.fillText(tagName(s.at),cx+Math.cos(a)*(L+15),cy-Math.sin(a)*(L+15));
  });
  // 前碳：实心点 + 三条线
  subOf(frontId).forEach(s=>{
    const a=Math.atan2(s.y,s.x);
    const x1=cx+Math.cos(a)*L,y1=cy-Math.sin(a)*L;
    ctx.strokeStyle=tagCol(s.at.tag);ctx.lineWidth=2.4;
    ctx.beginPath();ctx.moveTo(cx,cy);ctx.lineTo(x1,y1);ctx.stroke();
    ctx.fillStyle=tagCol(s.at.tag);ctx.font='600 13px -apple-system,sans-serif';
    ctx.textAlign='center';ctx.textBaseline='middle';
    ctx.fillText(tagName(s.at),cx+Math.cos(a)*(L+15),cy-Math.sin(a)*(L+15));
  });
  const g=ctx.createRadialGradient(cx,cy,1,cx,cy,gap-3);
  g.addColorStop(0,'#dfe9f7');g.addColorStop(1,'#4a5a76');
  ctx.fillStyle=g;ctx.beginPath();ctx.arc(cx,cy,gap-3,0,Math.PI*2);ctx.fill();
  ctx.strokeStyle='#0c121c';ctx.lineWidth=1.2;ctx.stroke();
}

/* ============================================================
   能量曲线图
============================================================ */
class EnergyChart{
  constructor(canvasId,opt){
    this.cv=document.getElementById(canvasId);
    this.opt=opt; this.x=opt.x0??0;
    this.cv.addEventListener('pointerdown',e=>{this._drag=true;this._pick(e);});
    window.addEventListener('pointermove',e=>{if(this._drag)this._pick(e);});
    window.addEventListener('pointerup',()=>this._drag=false);
    new ResizeObserver(()=>this.draw()).observe(this.cv);
  }
  _pick(e){
    const r=this.cv.getBoundingClientRect();
    const padL=38,padR=14,padT=22,padB=30;
    const u=(e.clientX-r.left-padL)/(r.width-padL-padR);
    const x=this.opt.xMin+clamp(u,0,1)*(this.opt.xMax-this.opt.xMin);
    this.opt.onSeek(x);
  }
  setX(x){this.x=x;this.draw();}
  draw(){
    const dpr=Math.min(devicePixelRatio,2);
    const W=this.cv.clientWidth,H=this.cv.clientHeight;
    if(!W)return;
    this.cv.width=W*dpr;this.cv.height=H*dpr;
    const ctx=this.cv.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);
    const padL=38,padR=14,padT=24,padB=30;
    const pw=W-padL-padR,ph=H-padT-padB;
    const {xMin,xMax,yMax}=this.opt;
    const X=x=>padL+(x-xMin)/(xMax-xMin)*pw;
    const Y=y=>padT+ph-y/yMax*ph;
    // 网格
    ctx.strokeStyle='rgba(120,150,190,.15)';ctx.fillStyle='#66799a';
    ctx.font='10.5px -apple-system,sans-serif';ctx.lineWidth=1;
    [0,.25,.5,.75,1].forEach(f=>{const y=padT+ph*f;ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(W-padR,y);ctx.stroke();
      ctx.textAlign='right';ctx.fillText(Math.round(yMax*(1-f)),padL-5,y+3);});
    // 曲线
    ctx.beginPath();
    for(let i=0;i<=240;i++){const x=xMin+(xMax-xMin)*i/240,y=this.opt.xToY(x);
      i?ctx.lineTo(X(x),Y(y)):ctx.moveTo(X(x),Y(y));}
    const grad=ctx.createLinearGradient(0,padT,0,padT+ph);
    grad.addColorStop(0,'#7fb4ff');grad.addColorStop(1,'#46d6c8');
    ctx.strokeStyle=grad;ctx.lineWidth=2.4;ctx.stroke();
    // 关键点
    (this.opt.keys||[]).forEach(k=>{
      const y=this.opt.xToY(k.x);
      ctx.beginPath();ctx.arc(X(k.x),Y(y),4.2,0,Math.PI*2);
      ctx.fillStyle=k.type==='bad'?'#ff6378':k.type==='mid'?'#a48bff':'#3ad79b';ctx.fill();
      ctx.font='10.5px -apple-system,"PingFang SC",sans-serif';ctx.textAlign='center';
      ctx.fillStyle='#9db2d2';
      const yy = k.labelBelow? Y(y)+15 : Y(y)-9;
      ctx.fillText(k.name,X(k.x),yy);
    });
    // x 轴标签
    ctx.fillStyle='#66799a';ctx.textAlign='center';
    (this.opt.xLabels||[]).forEach(l=>ctx.fillText(l.t,X(l.x),H-9));
    // 当前位置
    const cx=X(this.x),cy=Y(this.opt.xToY(this.x));
    ctx.setLineDash([4,4]);ctx.strokeStyle='rgba(220,232,255,.55)';ctx.lineWidth=1.2;
    ctx.beginPath();ctx.moveTo(cx,padT);ctx.lineTo(cx,padT+ph);ctx.stroke();ctx.setLineDash([]);
    ctx.beginPath();ctx.arc(cx,cy,5.5,0,Math.PI*2);ctx.fillStyle='#fff';ctx.fill();
    ctx.beginPath();ctx.arc(cx,cy,3.2,0,Math.PI*2);ctx.fillStyle='#2f7cf0';ctx.fill();
  }
}

/* ============================================================
   初始化各查看器
============================================================ */
const viewers={};
function initViewers(){
  viewers.ethane=new Viewer('vw-ethane',{r:6.6,phi:1.12,theta:.9,tx:0,ty:0,tz:.77});
  viewers.butane=new Viewer('vw-butane',{r:7.4,phi:1.12,theta:.9,tx:0,ty:0,tz:.77});
  viewers.rings=new Viewer('vw-rings',{r:6.2,phi:1.05,theta:.6,tx:0,ty:0,tz:.1});
  viewers.rings.opts.angleLabel=true;
  viewers.chair=new Viewer('vw-chair',{r:7.6,phi:1.18,theta:.7,tx:0,ty:0,tz:.15});
  viewers.axeq=new Viewer('vw-axeq',{r:7.8,phi:1.15,theta:.7,tx:0,ty:0,tz:.15});
  viewers.disub=new Viewer('vw-disub',{r:8.6,phi:1.15,theta:.7,tx:0,ty:0,tz:.15});
}

/* ============================================================
   工具栏通用
============================================================ */
function syncVdwBtn(btn,vw){
  const s=vw.opts.vdwScale;
  const label=s<.05?'关':s<.85?'半':'实';
  btn.textContent='◯ 填充:'+label;
  btn.classList.toggle('on',s>.05);
}
function bindToolbar(id,vw,onChange){
  document.getElementById(id).querySelectorAll('[data-toggle]').forEach(btn=>{
    btn.addEventListener('click',()=>{
      const key=btn.dataset.toggle;
      const map={vdw:'vdw',rot:'rotate',label:'labels',angle:'angleLabel',flag:'flags',aecolor:'aeColor'};
      const optKey=map[key];
      if(optKey==='vdw'){
        // 三档循环:关(0) → 半(.6) → 实(1)
        const cur=vw.opts.vdwScale;
        const nxt=cur<.05?.6:cur<.85?1:0;
        vw.opts.vdwScale=nxt;
        vw.opts.vdw=nxt>0;
        syncVdwBtn(btn,vw);
      }else{
        vw.opts[optKey]=!vw.opts[optKey];
        btn.classList.toggle('on',vw.opts[optKey]);
      }
      onChange&&onChange();
    });
    // 初始同步 vdw 按钮文字
    if(btn.dataset.toggle==='vdw') syncVdwBtn(btn,vw);
  });
  document.getElementById(id).querySelectorAll('[data-act]').forEach(btn=>{
    btn.addEventListener('click',()=>{
      if(btn.dataset.act==='resetView')vw.resetView();
      if(btn.dataset.act==='flip')onChange&&onChange('flip');
    });
  });
}

/* ============================================================
   模块一：乙烷
============================================================ */
const eth={angle:60,playing:false};
function ethaneEnergy(p){return 6.3*(1+Math.cos(3*p*D2R));}
function renderEthane(){
  const mol=buildEthane(eth.angle);
  viewers.ethane.setMolecule(mol);
  drawNewman(document.getElementById('nm-ethane'),mol,mol.front,mol.back);
  charts.ethane.setX(eth.angle);
  document.getElementById('eth-hh').textContent=(mol.info.hh*100).toFixed(0)+' pm';
  document.getElementById('eth-hh').style.color=mol.info.hh>=2.4?'#3ad79b':'#ff6378';
  const E=ethaneEnergy(eth.angle);
  document.getElementById('eth-energy').innerHTML=E.toFixed(1)+' <small>kJ/mol</small>';
  const ecl=angDist(eth.angle%120,0)<18;
  document.getElementById('eth-name').textContent=ecl?'重叠式':'交叉式';
  document.getElementById('eth-sname').textContent=ecl?'重叠式 Eclipsed':'交叉式 Staggered';
  document.getElementById('eth-desc').textContent=ecl?'前后 C–H 键两两正对，氢距最近':'前后氢彼此错开，相距最远';
  document.getElementById('eth-strains').innerHTML=ecl
    ?'<span class="strain-pill hot">H···H 距离 &lt; 240 pm</span><span class="strain-pill hot">扭转张力最大</span>'
    :'<span class="strain-pill cool">无空间张力</span><span class="strain-pill cool">扭转张力最小</span>';
  document.querySelectorAll('[data-eth]').forEach(b=>b.classList.toggle('on',angDist(eth.angle,+b.dataset.eth)<12));
  document.getElementById('eth-angle').value=eth.angle;
}
function initEthane(){
  bindToolbar('tb-ethane',viewers.ethane,renderEthane);
  const slider=document.getElementById('eth-angle');
  slider.addEventListener('input',()=>{stopPlay('eth');eth.angle=+slider.value;renderEthane();});
  document.querySelectorAll('[data-eth]').forEach(b=>b.addEventListener('click',()=>{
    stopPlay('eth');const target=+b.dataset.eth,start=eth.angle,delta=target-start;
    tweenTo(380,t=>{eth.angle=((start+delta*t)%360+360)%360;renderEthane();});
  }));
  document.getElementById('eth-play').addEventListener('click',function(){
    this.classList.toggle('on');togglePlay('eth',this,'连续旋转');
  });
}

/* ============================================================
   模块二：丁烷
============================================================ */
const but={angle:180,playing:false};
const BUT_KEYS=[[0,19],[60,3.8],[120,16],[180,0],[240,16],[300,3.8],[360,19]];
function butaneEnergy(t){return csLerp(BUT_KEYS,t);}
function butStateName(t){
  const cands=[[0,'全重叠式'],[60,'邻位交叉式'],[120,'部分重叠式'],[180,'对位交叉式'],[240,'部分重叠式'],[300,'邻位交叉式'],[360,'全重叠式']];
  let best=cands[0],bd=999;
  cands.forEach(c=>{const d=angDist(t,c[0]);if(d<bd){bd=d;best=c;}});
  return [best[1],bd];
}
function renderButane(){
  const mol=buildButane(but.angle);
  viewers.butane.setMolecule(mol);
  drawNewman(document.getElementById('nm-butane'),mol,mol.front,mol.back);
  charts.butane.setX(but.angle);
  const mm=mol.info.mm;
  document.getElementById('but-mm').innerHTML=(mm*100).toFixed(0)+' <small>pm</small>';
  document.getElementById('eth-mm').textContent=(mm*100).toFixed(0)+' pm';
  document.getElementById('eth-mm').style.color=mm>=3.4?'#3ad79b':mm>=2.9?'#c9a7ff':'#ff6378';
  const E=butaneEnergy(but.angle);
  document.getElementById('but-energy').innerHTML=E.toFixed(1)+' <small>kJ/mol</small>';
  const [name,d]=butStateName(but.angle);
  const en={全重叠式:'Fully eclipsed',邻位交叉式:'Gauche',部分重叠式:'Partially eclipsed',对位交叉式:'Anti'}[name];
  document.getElementById('but-name').textContent=name;
  document.getElementById('but-sname').textContent=name.replace('式','')+'式 '+en;
  const descs={全重叠式:'两个甲基正对，电子云强烈交叠',邻位交叉式:'两个甲基相邻，轻微空间位阻',
    部分重叠式:'甲基与氢重叠，斥力较大',对位交叉式:'两个甲基处于对位，相距最远'};
  document.getElementById('but-desc').textContent=descs[name];
  const pills={
    全重叠式:'<span class="strain-pill hot">甲基强烈空间位阻</span><span class="strain-pill hot">重叠式扭转张力</span>',
    部分重叠式:'<span class="strain-pill warm">Me···H 重叠斥力</span><span class="strain-pill warm">扭转张力较大</span>',
    邻位交叉式:'<span class="strain-pill warm">甲基轻微靠近</span><span class="strain-pill cool">交叉式</span>',
    对位交叉式:'<span class="strain-pill cool">空间位阻最小</span><span class="strain-pill cool">交叉式 · 优势构象</span>'};
  document.getElementById('but-strains').innerHTML=pills[name];
  document.querySelectorAll('[data-but]').forEach(b=>b.classList.toggle('on',angDist(but.angle,+b.dataset.but)<18));
  document.getElementById('but-angle').value=but.angle;
}
function initButane(){
  bindToolbar('tb-butane',viewers.butane,renderButane);
  const slider=document.getElementById('but-angle');
  slider.addEventListener('input',()=>{stopPlay('but');but.angle=+slider.value;renderButane();});
  document.querySelectorAll('[data-but]').forEach(b=>b.addEventListener('click',()=>{
    stopPlay('but');const target=+b.dataset.but,start=but.angle;
    let delta=target-start;
    tweenTo(380,t=>{but.angle=(start+delta*t);renderButane();});
  }));
  document.getElementById('but-play').addEventListener('click',function(){
    this.classList.toggle('on');togglePlay('but',this,'连续旋转');
  });
}

/* ============================================================
   模块三：环烷烃
============================================================ */
const ring={n:3};
const RING_INFO={
  3:{name:'环丙烷 Cyclopropane',desc:'三碳必须共面成等边三角形',strain:'很大',
    pills:'<span class="strain-pill hot">角张力大（60°）</span><span class="strain-pill hot">弯曲键（香蕉键）</span><span class="strain-pill warm">重叠式扭转张力</span>'},
  4:{name:'环丁烷 Cyclobutane',desc:'折叠成蝴蝶式以缓解扭转张力',strain:'较大',
    pills:'<span class="strain-pill hot">角张力较大（约 88°）</span><span class="strain-pill warm">仍有扭转张力</span>'},
  5:{name:'环戊烷 Cyclopentane',desc:'信封式：一碳翘出平面，相邻氢取交叉式',strain:'很小',
    pills:'<span class="strain-pill warm">轻微角张力（约 105°）</span><span class="strain-pill cool">扭转张力小</span>'},
  6:{name:'环己烷椅式 Cyclohexane (chair)',desc:'键角保持约 109.5°，全部 C–H 交叉式',strain:'无',
    pills:'<span class="strain-pill cool">无角张力</span><span class="strain-pill cool">无扭转张力</span><span class="strain-pill cool">无空间张力</span>'}
};
function renderRings(){
  const mol=buildRing(ring.n);
  viewers.rings.setMolecule(mol);
  const info=RING_INFO[ring.n];
  document.getElementById('ring-name').textContent=info.name.split(' ')[0];
  document.getElementById('ring-sname').textContent=info.name;
  document.getElementById('ring-desc').textContent=info.desc;
  document.getElementById('ring-strains').innerHTML=info.pills;
  const ang=mol.info.avgAngle;
  document.getElementById('ring-angleval').innerHTML=Math.round(ang)+'° <small style="color:var(--mut)">(理想 109.5°)</small>';
  document.getElementById('ring-angle-badge').innerHTML='平均键角：<b>'+Math.round(ang)+'°</b>';
  const se=document.getElementById('ring-strain');
  se.textContent=info.strain;
  se.style.color=ring.n===3?'var(--bad)':ring.n===4?'#ff9a7a':ring.n===5?'#a48bff':'var(--ok)';
  document.querySelectorAll('[data-ring]').forEach(b=>b.classList.toggle('on',+b.dataset.ring===ring.n));
}
function initRings(){
  bindToolbar('tb-rings',viewers.rings,renderRings);
  document.querySelectorAll('[data-ring]').forEach(b=>b.addEventListener('click',()=>{
    ring.n=+b.dataset.ring; renderRings();
  }));
}

/* ============================================================
   模块四：椅式 vs 船式
============================================================ */
const chair={t:0,playing:false,dir:1};
function chairEnergy(t){return csLerp(CHAIR_E_KEYS,t*100);}
function chairName(t){
  let best=CHAIR_NAME_KEYS[0],bd=999;
  CHAIR_NAME_KEYS.forEach(k=>{const d=Math.abs(k[0]-t*100);if(d<bd){bd=d;best=k;}});
  return best[1];
}
function renderChair(){
  const mol=buildChairMorph(chair.t);
  viewers.chair.setMolecule(mol);
  drawNewman(document.getElementById('nm-chair'),mol,mol.front,mol.back);
  charts.chair.setX(chair.t*100);
  const fl=mol.info.flag*100;
  document.getElementById('chair-flagval').innerHTML=mol.info.boatness>0.12?fl.toFixed(0)+' <small>pm</small>':'— <small>pm</small>';
  document.getElementById('chair-flag-badge').innerHTML=mol.info.boatness>0.12
    ?'旗杆 H···H：<b style="color:'+(mol.info.flag<2.4?'#ff6378':'#3ad79b')+'">'+fl.toFixed(0)+' pm</b>'
    :'旗杆 H···H：<b style="color:#8fa0b8">无（椅式氢全部错开）</b>';
  const E=chairEnergy(chair.t);
  document.getElementById('chair-energy').innerHTML=E.toFixed(1)+' <small>kJ/mol</small>';
  const nm=chairName(chair.t);
  const cn={椅式:'椅式 Chair','椅式′':'翻环椅式 Chair',半椅式:'半椅式 Half-chair',扭船式:'扭船式 Twist-boat',船式:'船式 Boat'}[nm];
  document.getElementById('chair-name').textContent=nm.replace('′','（翻转）');
  document.getElementById('chair-sname').textContent=cn;
  const nearBoat=chair.t>.25&&chair.t<.75;
  const atBoat=Math.abs(chair.t-.5)<.1, atTwist=Math.abs(chair.t-.34)<.08||Math.abs(chair.t-.66)<.08;
  document.getElementById('chair-desc').textContent=
    !nearBoat?'无角张力，无空间与扭转张力':
    atBoat?'旗杆氢强烈相斥，船底 C–H 全重叠':
    atTwist?'扭转使旗杆氢错开，斥力部分缓解':'翻环过程中的最高能量构象';
  let pills;
  if(!nearBoat) pills='<span class="strain-pill cool">无角张力</span><span class="strain-pill cool">无空间张力</span><span class="strain-pill cool">无扭转张力</span>';
  else if(atBoat) pills='<span class="strain-pill cool">无角张力</span><span class="strain-pill hot">旗杆氢空间张力（183 pm）</span><span class="strain-pill hot">重叠式扭转张力</span>';
  else if(atTwist) pills='<span class="strain-pill cool">无角张力</span><span class="strain-pill warm">空间张力缓解</span><span class="strain-pill warm">扭转张力缓解</span>';
  else pills='<span class="strain-pill hot">角/扭转张力最大（能垒）</span>';
  document.getElementById('chair-strains').innerHTML=pills;
  document.getElementById('chair-morph').value=chair.t*100;
  document.querySelectorAll('[data-morph]').forEach(b=>b.classList.toggle('on',Math.abs(chair.t*100-(+b.dataset.morph))<4));
}
function initChair(){
  bindToolbar('tb-chair',viewers.chair,renderChair);
  const slider=document.getElementById('chair-morph');
  slider.addEventListener('input',()=>{stopPlay('chair');chair.t=+slider.value/100;renderChair();});
  document.querySelectorAll('[data-morph]').forEach(b=>b.addEventListener('click',()=>{
    stopPlay('chair');const target=+b.dataset.morph/100,start=chair.t;
    tweenTo(520,tt=>{chair.t=lerp(start,target,tt);renderChair();});
  }));
  document.getElementById('chair-play').addEventListener('click',function(){
    this.classList.toggle('on');togglePlay('chair',this,'翻环动画');
  });
}

/* ============================================================
   模块五：a / e 键
============================================================ */
const axeq={sub:'Me',pos:'e',t:0,flip:0};
function axeqSubs(t){
  if(axeq.sub==='H')return [];
  // C1 = k0；e 起始为 up-equatorial? k0 偶数：a-up, e-down。演示“取代基朝上”：a=up, e 用 down 也可以，
  // 但为了一取代直接对照，令 a 位置取 up，e 位置取 down（视觉上都在 C1 外侧）
  return [{k:0,kind:axeq.sub,face:axeq.pos==='a'?'up':'down'}];
}
function renderAxeq(){
  const subs=axeqSubs(axeq.t);
  const mol=buildSubChair(subs,axeq.t);
  viewers.axeq.setMolecule(mol);
  // 统计取代基相关的近距接触
  const subIds=new Set();
  mol.atoms.forEach((a,i)=>{if(a._sub)subIds.add(i);});
  let nContacts=0;
  viewers.axeq.clashes.forEach(c=>{
    if(subIds.has(c.i)||subIds.has(c.j)) nContacts++;
  });
  const atEnd=axeq.t<.05||axeq.t>.95;
  const nowPos = axeq.t<.5 ? axeq.pos : (axeq.pos==='a'?'e':'a');
  const A=axeq.sub==='tBu'?21:7.6;
  const startAx=axeq.pos==='a'?1:0, endAx=axeq.pos==='a'?0:1;
  const axialness=lerp(startAx,endAx,smooth(axeq.t));
  const penalty=(axeq.sub==='H'?0:A)*axialness;
  document.getElementById('axeq-penalty').innerHTML=Math.round(penalty*10)/10+' <small>kJ/mol</small>';
  document.getElementById('axeq-n').innerHTML=nContacts+' <small>对</small>';
  document.getElementById('axeq-contact').innerHTML='取代基近距接触：<b>'+nContacts+'</b> 对';
  const subName={H:'氢（无取代）',Me:'甲基',tBu:'叔丁基'}[axeq.sub];
  document.getElementById('axeq-badge').textContent='取代基：'+subName+(atEnd?' · '+nowPos+'键':'（翻环中）');
  let title,desc,pills;
  if(axeq.sub==='H'){title='环己烷的 12 根 C–H 键';desc='6 根 a 键与 6 根 e 键，各 3 上 3 下';
    pills='<span class="strain-pill cool">a 键：竖直平行环轴</span><span class="strain-pill cool">e 键：斜伸环外</span>';}
  else if(!atEnd){title='翻环进行中…';desc='a、e 键正在互换，朝上/朝下取向不变';
    pills='<span class="strain-pill warm">a ⇄ e 互换中</span>';}
  else if(nowPos==='e'){title='e 键取代（优势构象）';desc='大基团伸向环外侧，无拥挤';
    pills='<span class="strain-pill cool">基团指向环外</span><span class="strain-pill cool">无 1,3-二直立斥力</span>';}
  else{title='a 键取代（非优势构象）';desc='与 C3、C5 的直立氢发生 1,3-二直立相互作用';
    pills='<span class="strain-pill hot">1,3-二直立空间位阻</span><span class="strain-pill warm">能量 +'+A+' kJ/mol</span>';}
  document.getElementById('axeq-sname').textContent=title;
  document.getElementById('axeq-desc').textContent=desc;
  document.getElementById('axeq-strains').innerHTML=pills;
}
function flipAxeq(){
  stopTweens();const start=axeq.t,end=axeq.t>.5?0:1;
  tweenTo(900,t=>{axeq.t=lerp(start,end,t);renderAxeq();},()=>{axeq.flip=end;});
}
function initAxeq(){
  bindToolbar('tb-axeq',viewers.axeq,act=>{if(act==='flip')flipAxeq();else renderAxeq();});
  document.querySelectorAll('#axeq-sub-seg button').forEach(b=>b.addEventListener('click',()=>{
    axeq.sub=b.dataset.sub;axeq.t=0;
    document.querySelectorAll('#axeq-sub-seg button').forEach(x=>x.classList.toggle('on',x===b));
    renderAxeq();
  }));
  document.querySelectorAll('#axeq-pos-seg button').forEach(b=>b.addEventListener('click',()=>{
    axeq.pos=b.dataset.pos;axeq.t=0;
    document.querySelectorAll('#axeq-pos-seg button').forEach(x=>x.classList.toggle('on',x===b));
    renderAxeq();
  }));
}

/* ============================================================
   模块六：二取代环己烷
============================================================ */
const disub={pos:'1,2',ct:'trans',t:0,tbu:false};
function disubSubs(){
  const j=+disub.pos.split(',')[1]-1;
  const face0='up';
  const face1=disub.ct==='cis'?'up':'down';
  return [{k:0,kind:'Me',face:face0},{k:j,kind:disub.tbu?'tBu':'Me',face:face1}];
}
function slotAt(k,face,flip){
  const f=analyticChair(flip).dirs[k];
  if(face==='up') return f.aUp?'a':'e';
  return f.aUp?'e':'a';
}
function renderDisub(){
  const subs=disubSubs();
  const mol=buildSubChair(subs,disub.t);
  viewers.disub.setMolecule(mol);
  const flip=disub.t>.5?1:0;
  const types=subs.map(s=>slotAt(s.k,s.face,flip));
  const combo=types.join('');
  const j=+disub.pos.split(',')[1];
  const cnName=(disub.ct==='cis'?'顺式-':'反式-')+disub.pos+'-二'+(disub.tbu?'甲基-4-叔丁基':'甲基')+'环己烷';
  document.getElementById('disub-sname').textContent=disub.tbu?'顺-1-甲基-4-叔丁基环己烷':cnName;
  document.getElementById('disub-badge').textContent=disub.tbu?'顺-1-甲基-4-叔丁基':'顺/反-'+disub.pos+'二甲基';
  document.getElementById('disub-combo-badge').innerHTML='当前构象：<b>'+combo+'</b>';
  document.getElementById('disub-combo').textContent='当前椅式：'+combo.split('').join('、')+' 键'+(flip?'（翻环后）':'');
  const nAxial=types.filter(x=>x==='a').length;
  const penalty=subs.reduce((sum,s)=>{const ty=slotAt(s.k,s.face,flip);return sum+(ty==='a'?(s.kind==='tBu'?21:7.6):0);},0);
  document.getElementById('disub-na').innerHTML=nAxial+' <small>个</small>';
  document.getElementById('disub-energy').innerHTML=penalty.toFixed(1)+' <small>kJ/mol</small>';
  let pills;
  if(disub.t<.05||disub.t>.95){
    if(nAxial===0)pills='<span class="strain-pill cool">ee · 双 e 键优势构象</span>';
    else if(nAxial===2)pills='<span class="strain-pill hot">aa · 双 a 键，翻环后得 ee</span>';
    else pills='<span class="strain-pill warm">ea/ae · 总有一个直立基团</span>';
    if(disub.tbu){
      const tbuType=slotAt(j-1,'up',flip);
      pills+=tbuType==='e'
        ?'<span class="strain-pill cool">叔丁基锁定 e 键</span>'
        :'<span class="strain-pill hot">叔丁基在 a 键（极少存在）</span>';
    }
  }else pills='<span class="strain-pill warm">翻环中：顺反不变，a ⇄ e 互换</span>';
  document.getElementById('disub-strains').innerHTML=pills;
}
function flipDisub(){
  stopTweens();const start=disub.t,end=disub.t>.5?0:1;
  tweenTo(900,t=>{disub.t=lerp(start,end,t);renderDisub();});
}
function initDisub(){
  bindToolbar('tb-disub',viewers.disub,act=>{if(act==='flip')flipDisub();});
  document.querySelectorAll('[data-pos12]').forEach(b=>b.addEventListener('click',()=>{
    disub.pos=b.dataset.pos12;disub.t=0;disub.tbu=false;
    document.querySelectorAll('[data-pos12]').forEach(x=>x.classList.toggle('on',x===b));
    renderDisub();
  }));
  document.querySelectorAll('#disub-ct-seg button').forEach(b=>b.addEventListener('click',()=>{
    disub.ct=b.dataset.ct;disub.t=0;disub.tbu=false;
    document.querySelectorAll('#disub-ct-seg button').forEach(x=>x.classList.toggle('on',x===b));
    renderDisub();
  }));
  document.getElementById('load-tbu').addEventListener('click',()=>{
    disub.tbu=true;disub.pos='1,4';disub.ct='cis';disub.t=0;
    document.querySelectorAll('[data-pos12]').forEach(x=>x.classList.toggle('on',x.dataset.pos12==='1,4'));
    document.querySelectorAll('#disub-ct-seg button').forEach(x=>x.classList.toggle('on',x.dataset.ct==='cis'));
    renderDisub();
  });
}

/* ============================================================
   连续旋转播放
============================================================ */
const plays={};
function togglePlay(name,btn,label){
  plays[name]=!plays[name];
  if(plays[name]){
    btn.textContent='⏸ 暂停';
    btn._last=performance.now();
    ensurePump();
  }else{btn.textContent='▶ '+label;}
}
function stopPlay(name){
  plays[name]=false;
  const map={eth:['eth-play','连续旋转'],but:['but-play','连续旋转'],chair:['chair-play','翻环动画']};
  if(map[name]){const b=document.getElementById(map[name][0]);b.classList.remove('on');b.textContent='▶ '+map[name][1];}
}
function tickPlays(now){
  const dt=Math.min(50,now-(plays._last??now));plays._last=now;
  if(plays.eth){eth.angle=(eth.angle+dt*.06)%360;renderEthane();}
  if(plays.but){but.angle=(but.angle+dt*.06)%360;renderButane();}
  if(plays.chair){
    chair.t+=dt*.00028*chair.dir;
    if(chair.t>=1){chair.t=1;chair.dir=-1;}
    if(chair.t<=0){chair.t=0;chair.dir=1;}
    renderChair();
  }
}

/* ============================================================
   导航与进度保存
============================================================ */
const STORE_KEY='conformation-lab-v1';
const store=JSON.parse(localStorage.getItem(STORE_KEY)||'{}');
store.understood=store.understood||{};
const PAGES=['intro','ethane','butane','rings','chair','axeq','disub'];
const PAGE_RENDER={ethane:renderEthane,butane:renderButane,rings:renderRings,chair:renderChair,axeq:renderAxeq,disub:renderDisub};
function showPage(name){
  document.querySelectorAll('.page').forEach(p=>p.classList.toggle('show',p.id==='page-'+name));
  document.querySelectorAll('.nav-item').forEach(b=>b.classList.toggle('active',b.dataset.page===name));
  store.lastPage=name;localStorage.setItem(STORE_KEY,JSON.stringify(store));
  document.getElementById('main').scrollTop=0;
  // display 切换后重新测量尺寸并重绘
  requestAnimationFrame(()=>{
    Object.values(viewers).forEach(v=>v._resize());
    if(PAGE_RENDER[name])PAGE_RENDER[name]();
  });
}
function updateProgress(){
  const n=PAGES.filter(p=>store.understood[p]).length;
  document.getElementById('progBar').style.width=(n/PAGES.length*100)+'%';
  document.getElementById('progPct').textContent=Math.round(n/PAGES.length*100);
  PAGES.forEach(p=>document.querySelector(`.nav-item[data-page="${p}"]`)?.classList.toggle('done',!!store.understood[p]));
  document.querySelectorAll('[data-understand]').forEach(b=>{
    const on=!!store.understood[b.dataset.understand];
    b.classList.toggle('done',on);
    b.textContent=on?'✓ 已理解':'我已理解本节';
  });
}
function initNav(){
  document.querySelectorAll('.nav-item').forEach(b=>b.addEventListener('click',()=>showPage(b.dataset.page)));
  document.querySelectorAll('[data-understand]').forEach(b=>b.addEventListener('click',()=>{
    const k=b.dataset.understand;
    store.understood[k]=!store.understood[k];
    localStorage.setItem(STORE_KEY,JSON.stringify(store));
    updateProgress();
  }));
}

/* ============================================================
   启动
============================================================ */
const charts={};
function initCharts(){
  charts.ethane=new EnergyChart('ch-ethane',{
    xMin:0,xMax:360,yMax:14,x0:60,xToY:ethaneEnergy,onSeek:x=>{stopPlay('eth');eth.angle=x;renderEthane();},
    keys:[{x:0,name:'重叠',type:'bad'},{x:60,name:'交叉',type:'ok'},{x:120,name:'重叠',type:'bad'},{x:180,name:'交叉',type:'ok'},{x:240,name:'重叠',type:'bad'},{x:300,name:'交叉',type:'ok'}],
    xLabels:[{x:0,t:'0°'},{x:180,t:'180°'},{x:360,t:'360°'}]
  });
  charts.butane=new EnergyChart('ch-butane',{
    xMin:0,xMax:360,yMax:22,x0:180,xToY:butaneEnergy,onSeek:x=>{stopPlay('but');but.angle=x;renderButane();},
    keys:[{x:0,name:'全重叠',type:'bad'},{x:60,name:'邻位交叉',type:'mid'},{x:120,name:'部分重叠',type:'bad'},{x:180,name:'对位交叉',type:'ok',labelBelow:true},{x:240,name:'部分重叠',type:'bad'},{x:300,name:'邻位交叉',type:'mid'},{x:360,name:'全重叠',type:'bad',labelBelow:true}],
    xLabels:[{x:0,t:'0°'},{x:180,t:'180°'},{x:360,t:'360°'}]
  });
  charts.chair=new EnergyChart('ch-chair',{
    xMin:0,xMax:100,yMax:56,x0:0,xToY:x=>csLerp(CHAIR_E_KEYS,x),onSeek:x=>{stopPlay('chair');chair.t=x/100;renderChair();},
    keys:[{x:0,name:'椅式',type:'ok'},{x:17,name:'半椅',type:'bad'},{x:34,name:'扭船',type:'mid'},{x:50,name:'船式',type:'bad'},{x:66,name:'扭船',type:'mid'},{x:83,name:'半椅',type:'bad'},{x:100,name:'椅式′',type:'ok'}],
    xLabels:[{x:0,t:'椅式'},{x:50,t:'船式'},{x:100,t:'翻环椅式'}]
  });
}

initViewers();
initCharts();
initNav();
initEthane();initButane();initRings();initChair();initAxeq();initDisub();
renderEthane();renderButane();renderRings();renderChair();renderAxeq();renderDisub();
updateProgress();
showPage(PAGES.includes(store.lastPage)?store.lastPage:'intro');

/* ============================================================
   移动端布局重排:把 Newman/能量图 subviews 移到 grid 末尾,
   让 viewer → 工具栏 → 控制台(滑块+chip) 紧凑可见,提升互动性
============================================================ */
function applyMobileLayout(){
  const isMobile = window.innerWidth <= 768;
  document.querySelectorAll('.g-2').forEach(grid=>{
    const viewerCard = grid.children[0];
    if(!viewerCard) return;
    const subviews = viewerCard.querySelector(':scope > .subviews');
    if(!subviews) return;
    if(isMobile){
      if(subviews.parentElement === viewerCard) grid.appendChild(subviews);
    }else{
      if(subviews.parentElement !== viewerCard) viewerCard.appendChild(subviews);
    }
  });
  // 切换布局后 viewer 尺寸可能变化,触发重测
  Object.values(viewers).forEach(v=>v._resize && v._resize());
}
let _mobilRzT;
window.addEventListener('resize', ()=>{
  clearTimeout(_mobilRzT);
  _mobilRzT=setTimeout(applyMobileLayout, 150);
});
applyMobileLayout();

let lastTick=0;
function anyPlayActive(){return ['eth','but','chair'].some(k=>plays[k]);}
function loop(now){
  if(now-lastTick<16)return;              // 防止 rAF 与兜底定时器重复驱动
  lastTick=now;
  tickTweens(now);
  tickPlays(now);
  Object.values(viewers).forEach(v=>v.render());
  // 无活动动画时关闭兜底定时器，避免空闲期持续渲染
  if(pumpTimer!=null&&tweens.length===0&&!anyPlayActive()){
    clearInterval(pumpTimer);pumpTimer=null;
  }
}
requestAnimationFrame(function frame(now){loop(now);requestAnimationFrame(frame);});
