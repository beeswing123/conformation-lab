import * as THREE from './three.module.js';

/* ============================================================
   工具函数
============================================================ */
const V = (x,y,z)=>new THREE.Vector3(x,y,z);
const D2R = Math.PI/180;
const S3 = 2*Math.SQRT2/3;          // 0.9428 四面体方向水平分量
const CC = 1.54, CH = 1.09;         // 键长 Å
const CX = { F:1.39, Cl:1.78, Br:1.93, I:2.14, O:1.43 };  // C–杂原子键长
const VDW = { C:1.70, H:1.20, F:1.47, Cl:1.75, Br:1.85, I:1.98, O:1.52 };
// A 值:取代基处于 a 键的 1,3-二直立惩罚 (kJ/mol)，越大越偏好 e 键
const AVAL = { H:0, F:0.5, Cl:1.7, Br:1.8, I:1.9, OH:3.9, Me:7.6, Et:7.9, iPr:9.2, tBu:21 };
const GROUP_LABEL = { F:'F', Cl:'Cl', Br:'Br', OH:'OH', Me:'Me', Et:'Et', iPr:'i-Pr', tBu:'t-Bu' };
const GROUP_CN = { H:'氢', F:'氟', Cl:'氯', Br:'溴', OH:'羟基', Me:'甲基', Et:'乙基', iPr:'异丙基', tBu:'叔丁基' };
const GROUP_EN = { H:'H', F:'F', Cl:'Cl', Br:'Br', OH:'OH', Me:'Me', Et:'Et', iPr:'i-Pr', tBu:'t-Bu' };
let LANG='cn';   // 界面语言,需在所有渲染函数之前声明(函数内会读取)
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
  // 每个碳标一个代表性 H–C–H 键角(四面体 ~109.5°)
  const ang=(c,h1,h2)=>atoms[h1].p.clone().sub(atoms[c].p).angleTo(atoms[h2].p.clone().sub(atoms[c].p))/D2R;
  const angleData=[{at:0,n1:2,n2:3,v:ang(0,2,3)},{at:1,n1:5,n2:6,v:ang(1,5,6)}];
  return {atoms,bonds,measures,info:{hh},angleData, front:0, back:1};
}

// ---- 甲基（在中心碳 M 上，键轴 u 指向所连中心碳；tag 决定配色）----
function attachMethyl(atoms, bonds, M, u, tag='me', label){
  const e1 = Math.abs(u.y)>0.9 ? V(1,0,0) : V(0,1,0).sub(u.clone().multiplyScalar(u.y)).normalize();
  const e2 = new THREE.Vector3().crossVectors(u,e1).normalize();
  const mc = atoms.length;
  atoms.push({el:'C',p:M,tag,label:label||(tag==='tb'?'CH₃':'CH₃'),big:true});
  for(let i=0;i<3;i++){
    const a=i*120*D2R;
    const dir=u.clone().multiplyScalar(-1/3)
      .add(e1.clone().multiplyScalar(S3*Math.cos(a)))
      .add(e2.clone().multiplyScalar(S3*Math.sin(a))).normalize();
    atoms.push({el:'H',p:M.clone().add(dir.multiplyScalar(CH)),tag:tag+'H'});
    bonds.push([mc,atoms.length-1]);
  }
  return mc;
}
// ---- 叔丁基 ----
function attachTBu(atoms,bonds,S,u){
  const e1 = Math.abs(u.y)>0.9 ? V(1,0,0) : V(0,1,0).sub(u.clone().multiplyScalar(u.y)).normalize();
  const e2 = new THREE.Vector3().crossVectors(u,e1).normalize();
  const cen=atoms.length;
  atoms.push({el:'C',p:S,tag:'tb',label:'t-Bu',big:true});
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
// 通用 sp³ 碳:在 pos 放一个碳并补 n 个 H(n=1/2/3),键轴 u 指向所连环
const _basis=u=>{
  const e1=Math.abs(u.y)>.9?V(1,0,0):V(0,1,0).sub(u.clone().multiplyScalar(u.y)).normalize();
  return [e1,new THREE.Vector3().crossVectors(u,e1).normalize()];
};
function addCHn(atoms,bonds,pos,u,n,tag,lbl,angOff=0){
  const [e1,e2]=_basis(u);
  const ci=atoms.length;
  atoms.push({el:'C',p:pos,tag,label:lbl,big:true});
  for(let i=0;i<n;i++){
    const a=(angOff+i*(360/n))*D2R;
    const dir=u.clone().multiplyScalar(-1/3)
      .add(e1.clone().multiplyScalar(S3*Math.cos(a)))
      .add(e2.clone().multiplyScalar(S3*Math.sin(a))).normalize();
    atoms.push({el:'H',p:pos.clone().add(dir.multiplyScalar(CH)),tag:tag+'H'});
    bonds.push([ci,atoms.length-1]);
  }
  return ci;
}
// 环上取代基:k=环碳 idx,p=环碳位置,d=取代键单位方向(由环向外),kind=基团名
function attachGroup(atoms,bonds,k,p,d,kind){
  // 单原子卤素
  if(CX[kind]&&kind!=='O'){
    const ai=atoms.length;
    atoms.push({el:kind,p:p.clone().add(d.clone().multiplyScalar(CX[kind])),
      tag:kind,label:GROUP_LABEL[kind],big:true});
    bonds.push([k,ai]); return ai;
  }
  if(kind==='OH'){
    const [e1]=_basis(d);
    const oi=atoms.length;
    atoms.push({el:'O',p:p.clone().add(d.clone().multiplyScalar(CX.O)),tag:'OH',label:'OH',big:true});
    bonds.push([k,oi]);
    const hd=d.clone().multiplyScalar(-1/3).add(e1.clone().multiplyScalar(S3)).normalize();
    atoms.push({el:'H',p:atoms[oi].p.clone().add(hd.multiplyScalar(.96)),tag:'OHH'});
    bonds.push([oi,atoms.length-1]);
    return oi;
  }
  const u=d.clone().negate();
  const C1=p.clone().add(d.clone().multiplyScalar(CC));
  if(kind==='Me'){ const mc=attachMethyl(atoms,bonds,C1,u,'me'); atoms[mc].big=true; bonds.push([k,mc]); return mc; }
  if(kind==='tBu'){ const tc=attachTBu(atoms,bonds,C1,u); bonds.push([k,tc]); return tc; }
  if(kind==='Et'){
    const c1=addCHn(atoms,bonds,C1,u,2,'et','CH₂',30);
    bonds.push([k,c1]);
    const C2=C1.clone().add(d.clone().multiplyScalar(CC));
    const c2=attachMethyl(atoms,bonds,C2,d.clone().negate(),'et'); atoms[c2].big=true;
    bonds.push([c1,c2]); return c1;
  }
  if(kind==='iPr'){
    // 中心 CH:1 个 H + 2 个甲基,三等分绕 u
    const [e1,e2]=_basis(u);
    const c1=atoms.length;
    atoms.push({el:'C',p:C1,tag:'ip',label:'CH',big:true});
    bonds.push([k,c1]);
    const hdir=u.clone().multiplyScalar(-1/3).add(e1.clone().multiplyScalar(S3)).normalize();
    atoms.push({el:'H',p:C1.clone().add(hdir.multiplyScalar(CH)),tag:'ipH'});
    bonds.push([c1,atoms.length-1]);
    for(let i=0;i<2;i++){
      const a=((i?240:120))*D2R;
      const md=u.clone().multiplyScalar(-1/3).add(e1.clone().multiplyScalar(S3*Math.cos(a)))
        .add(e2.clone().multiplyScalar(S3*Math.sin(a))).normalize();
      const Mp=C1.clone().add(md.clone().multiplyScalar(CC));
      const mc=attachMethyl(atoms,bonds,Mp,md,'ip'); atoms[mc].big=true;
      bonds.push([c1,mc]);
    }
    return c1;
  }
  // 兜底:甲基
  const mc=attachMethyl(atoms,bonds,C1,u,'me'); bonds.push([k,mc]); return mc;
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
  return {atoms,bonds,measures,info:{mm},angleData:computeAngles(atoms,bonds),front:0,back:1};
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

// 通用键角:为有 ≥2 个非氢邻居的原子,取前两个(环/主链键在 bonds 中先入列)算 X–C–Y
function computeAngles(atoms,bonds){
  const adj=atoms.map(()=>[]);
  bonds.forEach(([a,b])=>{adj[a].push(b);adj[b].push(a);});
  const out=[];
  atoms.forEach((at,c)=>{
    if(at.el==='H')return;
    const heavy=adj[c].filter(n=>atoms[n].el!=='H');
    if(heavy.length<2)return;
    const n1=heavy[0],n2=heavy[1];
    const v=atoms[n1].p.clone().sub(at.p).angleTo(atoms[n2].p.clone().sub(at.p))/D2R;
    out.push({at:c,n1,n2,v});
  });
  return out;
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
  return {atoms,bonds,measures,info:{flag:best,boatness},angleData:computeAngles(atoms,bonds),front:1,back:2};
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
      const gi=attachGroup(atoms,bonds,k,p,d,s.kind);
      atoms[gi].slot=typeNow; atoms[gi]._sub=true;
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
  return {atoms,bonds,measures:[],info:{},angleData:computeAngles(atoms,bonds)};
}

/* ============================================================
   Three.js 分子查看器
============================================================ */
const COLORS={
  C:0xaeb9c9, H:0xeef3fb, me:0x5aa7ff, meH:0xb9dcff, tb:0x3fd4c0, tbH:0xbef5ee,
  et:0x46c7d6, etH:0xb6ecf2, ip:0xa48bff, ipH:0xd6c9ff,
  F:0x7be382, Cl:0x43c96b, Br:0xd9703a, I:0xb06bd9, OH:0xff6b78, OHH:0xeef3fb,
  ax:0x6ea8ff, eq:0x46d6c8, flag:0xff8fa0
};
// 元素文字标签颜色
const LBL_COLOR={me:'#8fc2ff',tb:'#7feadd',et:'#86ecf5',ip:'#c9b8ff',
  F:'#a6f0ad',Cl:'#86eca6',Br:'#ffb284',I:'#d3a8f5',OH:'#ff9aa5'};
const _spR={C:.31,H:.19,F:.24,Cl:.30,Br:.34,I:.40,O:.25};
const sphereGeo={},vdwGeo={};
Object.keys(VDW).forEach(el=>{
  sphereGeo[el]=new THREE.SphereGeometry(_spR[el]??.28,28,22);
  vdwGeo[el]=new THREE.SphereGeometry(VDW[el],26,20);
});
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
    this.scene.background=new THREE.Color('#0d141f');
    this._updateBg();   // 读取 CSS 变量 --vw-bg 设置场景背景
    this.camera=new THREE.PerspectiveCamera(42,1,0.1,100);
    this.cam0=cam0||{r:7.2,phi:1.18,theta:0.85,tx:0,ty:0,tz:0.1};
    this.cam={...this.cam0};
    this.opts={vdw:true,vdwScale:.6,rotate:true,labels:true,angleLabel:false,flags:false,aeColor:true,clash:true};
    this.panMode=false;   // 平移模式:开启后单指/左键拖动 = 平移
    this._cores=[];this._vdw=[];this._sprites=[];
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
    tweenTo(420,t=>{b.r=lerp(b.r,a.r,t);b.phi=lerp(b.phi,a.phi,t);b.theta=lerp(b.theta,a.theta,t);
      b.tx=lerp(b.tx,a.tx,t);b.ty=lerp(b.ty,a.ty,t);b.tz=lerp(b.tz,a.tz,t);}); }
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
  _updateBg(){
    // 从 CSS 变量 --vw-bg 读取背景色,日/夜切换时同步更新
    const v=getComputedStyle(document.body).getPropertyValue('--vw-bg').trim()||'#0d141f';
    this.scene.background=new THREE.Color(v);
  }
  _bind(){
    let drag=false,panning=false,lx=0,ly=0,pinch=0,mx=0,my=0;
    // 用户拖拽视角时,关闭"视角自转"按钮(让用户感知自己接管了)
    const takeOver=()=>{
      if(this.opts.rotate){
        this.opts.rotate=false;
        const card=this.canvas.closest('.card');
        const rotBtn=card?.querySelector('[data-toggle="rot"]');
        if(rotBtn) rotBtn.classList.remove('on');
      }
    };
    // 平移:模型跟手。dx>0(手指右)→ 视点左移;dy>0(手指下)→ 视点上移
    const pan=(dx,dy)=>{
      this.camera.updateMatrixWorld();
      const right=new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld,0);
      const up=new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld,1);
      const k=this.cam.r*.0018;
      this.cam.tx=clamp(this.cam.tx-right.x*dx*k+up.x*dy*k,-5,5);
      this.cam.ty=clamp(this.cam.ty-right.y*dx*k+up.y*dy*k,-4,4);
      this.cam.tz=clamp(this.cam.tz-right.z*dx*k+up.z*dy*k,-5,5);
    };
    const down=e=>{
      drag=true;
      panning=this.panMode||e.button===2;   // 平移模式 或 鼠标右键
      lx=e.clientX;ly=e.clientY;this._lastInteract=performance.now();
    };
    this.canvas.addEventListener('pointerdown',e=>{this.canvas.setPointerCapture(e.pointerId);down(e);});
    this.canvas.addEventListener('pointermove',e=>{
      if(!drag)return;
      const dx=e.clientX-lx,dy=e.clientY-ly;lx=e.clientX;ly=e.clientY;
      if(e.pointerType==='touch'&&this._d2)return;
      if(dx||dy) takeOver();
      if(panning){ pan(dx,dy); }
      else { this.cam.theta-=dx*.006; this.cam.phi=clamp(this.cam.phi-dy*.006,.15,Math.PI-.15); }
      this._lastInteract=performance.now();
    });
    const up=()=>{if(drag){this._lastInteract=performance.now();}drag=false;panning=false;};
    window.addEventListener('pointerup',up);
    // 右键拖动平移:阻止右键菜单
    this.canvas.addEventListener('contextmenu',e=>e.preventDefault());
    // 双击复位视角
    this.canvas.addEventListener('dblclick',()=>this.resetView());
    this.canvas.addEventListener('wheel',e=>{e.preventDefault();
      this.cam.r=clamp(this.cam.r*(1+e.deltaY*.0011),3.2,20);this._lastInteract=performance.now();},{passive:false});
    this.canvas.addEventListener('touchstart',e=>{
      if(e.touches.length===2){
        const dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;
        pinch=Math.hypot(dx,dy);
        mx=(e.touches[0].clientX+e.touches[1].clientX)/2;
        my=(e.touches[0].clientY+e.touches[1].clientY)/2;
        this._d2=true;drag=false;this._lastInteract=performance.now();
      }
    },{passive:true});
    this.canvas.addEventListener('touchmove',e=>{
      if(e.touches.length===2){
        const dx=e.touches[0].clientX-e.touches[1].clientX,dy=e.touches[0].clientY-e.touches[1].clientY;
        // 捏合缩放
        const d=Math.hypot(dx,dy);this.cam.r=clamp(this.cam.r*pinch/d,3.2,20);pinch=d;
        // 双指中点拖动 = 平移
        const cx=(e.touches[0].clientX+e.touches[1].clientX)/2;
        const cy=(e.touches[0].clientY+e.touches[1].clientY)/2;
        pan(cx-mx,cy-my); mx=cx; my=cy;
        this._lastInteract=performance.now();
      }
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
    // 原子:核心球 + VDW 球始终创建,显隐/配色交给 applyStyle 热更新(切换无卡顿)
    this._cores=[];this._vdw=[];this._sprites=[];
    mol.atoms.forEach((at,i)=>{
      const m=new THREE.Mesh(sphereGeo[at.el],new THREE.MeshPhongMaterial({shininess:55,specular:0x33445c}));
      m.position.copy(at.p); this.group.add(m); this._cores.push(m);
      const vm=new THREE.Mesh(vdwGeo[at.el],new THREE.MeshPhongMaterial({
        transparent:true,depthWrite:false,shininess:90}));
      vm.position.copy(at.p); this.group.add(vm); this._vdw.push(vm);
      if(at.label){
        const sp=makeSprite(at.label, LBL_COLOR[at.tag]||'#dce8fa');
        sp.position.copy(at.p).add(V(0,.58,0)); sp.scale.set(.96,.48,1);
        this.lg.add(sp); this._sprites.push({sp,kind:'label'});
      }
      if(at.tag==='ax'||at.tag==='eq'){
        const sp=makeSprite(at.tag==='ax'?'a':'e', at.tag==='ax'?'#8fbcff':'#76e7d6');
        sp.position.copy(at.p).add(V(0,.46,0)); sp.scale.set(.52,.42,1);
        this.lg.add(sp); this._sprites.push({sp,kind:'ae'});
      }
      if(at.tag==='flag'){
        const sp=makeSprite('旗杆H','#ffb3bd');
        sp.position.copy(at.p).add(V(0,.5,0)); sp.scale.set(.96,.46,1);
        this.lg.add(sp); this._sprites.push({sp,kind:'flag'});
      }
    });
    this.applyStyle();
    this.render();   // 立即上屏一帧（后台 rAF 暂停时也保证交互即时可见）
  }
  // 仅刷新显示样式(填充/标签/a-e/旗杆/红线),不重建分子 → 即时响应
  applyStyle(){
    const mol=this.mol; if(!mol)return;
    const s=this.opts;
    const tagOf=at=>{let t=at.tag;return (!s.aeColor&&(t==='ax'||t==='eq'))?'H':t;};
    this._cores.forEach((m,i)=>{
      const t=tagOf(mol.atoms[i]);
      m.material.color.setHex(COLORS[t]??COLORS[mol.atoms[i].el]);
      m.material.emissive.setHex(t==='flag'&&s.flags?0x55101b:0x000000);
    });
    this._vdw.forEach((vm,i)=>{
      const at=mol.atoms[i],t=tagOf(at);
      const on=s.vdw&&s.vdwScale>0; vm.visible=on;
      if(!on)return;
      const sev=this.atomSev.get(i)||0;
      const strong=s.clash&&(t==='flag'||sev>=.45);
      vm.material.color.setHex(strong?0xff4d66:(COLORS[t]??COLORS[at.el]));
      const opScale=Math.min(1,s.vdwScale+.18);
      vm.material.opacity=(strong?.22:(at.big?.16:.10))*opScale;
      vm.scale.setScalar(s.vdwScale===1?1:s.vdwScale);
    });
    this._sprites.forEach(o=>{
      o.sp.visible=o.kind==='label'?s.labels:o.kind==='ae'?s.aeColor:s.flags;
    });
    this.render();
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
    if(this.opts.clash) this.clashes.forEach(c=>{
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
  ctx.save();ctx.font='600 13.5px -apple-system,sans-serif';
  const w=ctx.measureText(text).width+16;
  ctx.fillStyle='rgba(13,20,32,.85)';ctx.strokeStyle=color;ctx.lineWidth=1.2;
  roundRect(ctx,x-w/2,y-11,w,22,7);ctx.fill();ctx.stroke();
  ctx.fillStyle=color;ctx.textAlign='center';ctx.textBaseline='middle';
  ctx.fillText(text,x,y+1);ctx.restore();
}
function roundRect(ctx,x,y,w,h,r){ctx.beginPath();ctx.moveTo(x+r,y);ctx.arcTo(x+w,y,x+w,y+h,r);
  ctx.arcTo(x+w,y+h,x,y+h,r);ctx.arcTo(x,y+h,x,y,r);ctx.arcTo(x,y,x+w,y,r);ctx.closePath();}

const spriteCache=new Map();
function makeSprite(text,color){
  const key=text+color;
  if(spriteCache.has(key))return spriteCache.get(key).clone();
  const c=document.createElement('canvas');c.width=320;c.height=160;
  const x=c.getContext('2d');
  x.font='700 66px -apple-system,"PingFang SC",sans-serif';x.textAlign='center';x.textBaseline='middle';
  x.shadowColor='rgba(0,0,0,.85)';x.shadowBlur=10;
  x.fillStyle=color;x.fillText(text,160,84);
  const tex=new THREE.CanvasTexture(c);tex.anisotropy=8;
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
  const en=typeof LANG!=='undefined'&&LANG==='en';
  const label=s<.05?(en?'Off':'关'):s<.85?(en?'Half':'半'):(en?'Full':'实');
  btn.textContent=(en?'◯ Fill:':'◯ 填充:')+label;
  btn.title=(en?'Space fill: ':'空间填充:')+(s<.05?(en?'Off':'关'):s<.85?(en?'Half':'半'):(en?'Full':'实'));
  btn.classList.toggle('on',s>.05);
}
// 纯显示类开关:只热更新样式,不重建分子(消除填充/标签切换卡顿)
const DISPLAY_TOGGLES={label:'labels',angle:'angleLabel',flag:'flags',aecolor:'aeColor',clash:'clash',rot:'rotate'};
// 超窄屏(≤480)工具栏改用纯图标单行显示,data-micon 供 CSS ::before 读取
const TB_ICON={vdw:'◉',rot:'⟳',pan:'✥',angle:'∠',clash:'⚠',label:'🏷',flag:'⚑',aecolor:'🎨'};
function bindToolbar(id,vw,onChange){
  document.getElementById(id).querySelectorAll('[data-toggle]').forEach(btn=>{
    const key=btn.dataset.toggle;
    btn.dataset.micon=TB_ICON[key]||'•';
    if(!btn.title) btn.title=btn.textContent.trim();
    if(key==='vdw') syncVdwBtn(btn,vw);
    else if(DISPLAY_TOGGLES[key]) btn.classList.toggle('on',!!vw.opts[DISPLAY_TOGGLES[key]]);
    btn.addEventListener('click',()=>{
      if(key==='vdw'){
        // 三档循环:关(0) → 半(.6) → 实(1)
        const cur=vw.opts.vdwScale;
        const nxt=cur<.05?.6:cur<.85?1:0;
        vw.opts.vdwScale=nxt; vw.opts.vdw=nxt>0;
        syncVdwBtn(btn,vw); vw.applyStyle();           // 即时,不重建
      }else if(key==='pan'){
        vw.panMode=!vw.panMode;
        btn.classList.toggle('on',vw.panMode);
        btn.textContent=vw.panMode?'✥ 平移中':'✥ 平移';
        btn.title=vw.panMode?'退出平移':'平移视图';
      }else if(DISPLAY_TOGGLES[key]){
        const ok=DISPLAY_TOGGLES[key];
        vw.opts[ok]=!vw.opts[ok];
        btn.classList.toggle('on',vw.opts[ok]);
        vw.applyStyle();                                // 键角/红线/标签即时刷新
      }
    });
  });
  document.getElementById(id).querySelectorAll('[data-act]').forEach(btn=>{
    btn.dataset.micon=btn.dataset.act==='flip'?'🔄':'↺';
    if(!btn.title) btn.title=btn.dataset.act==='flip'?'翻环':'复位视角';
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
  const A=AVAL[axeq.sub]||0;
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
const disub={pos:'1,2',ct:'trans',t:0,g1:'Me',g2:'Me'};
// 可选基团(按 A 值升序)
const DISUB_GROUPS=['F','Cl','Br','OH','Me','Et','iPr','tBu'];
function disubSubs(){
  const j=+disub.pos.split(',')[1]-1;
  const face1=disub.ct==='cis'?'up':'down';
  return [{k:0,kind:disub.g1,face:'up'},{k:j,kind:disub.g2,face:face1}];
}
function slotAt(k,face,flip){
  const f=analyticChair(flip).dirs[k];
  if(face==='up') return f.aUp?'a':'e';
  return f.aUp?'e':'a';
}
function disubName(){
  const cis=disub.ct==='cis';
  const n2=disub.pos.split(',')[1];
  const g1=GROUP_CN[disub.g1], g2=GROUP_CN[disub.g2];
  if(disub.g1===disub.g2) return `${cis?'顺':'反'}-${disub.pos}-${disub.g1==='Me'?'二甲基':'二'+g1}环己烷`;
  return `${cis?'顺':'反'}-1-${g1}-${n2}-${g2}环己烷`;
}
function renderDisub(){
  const subs=disubSubs();
  const mol=buildSubChair(subs,disub.t);
  viewers.disub.setMolecule(mol);
  const flip=disub.t>.5?1:0;
  const types=subs.map(s=>slotAt(s.k,s.face,flip));
  const combo=types.join('');
  const n2=disub.pos.split(',')[1];
  const title=document.getElementById('g2-title');
  if(title) title.innerHTML=`C${n2} 取代基 <small>（A 值 kJ/mol）</small>`;
  document.getElementById('disub-sname').textContent=disubName();
  document.getElementById('disub-badge').textContent=(disub.ct==='cis'?'顺':'反')+'-'+disub.pos+' '+GROUP_LABEL[disub.g1]+'/'+GROUP_LABEL[disub.g2];
  document.getElementById('disub-combo-badge').innerHTML='构象：<b>'+combo+'</b>';
  document.getElementById('disub-combo').textContent='当前椅式：'+combo.split('').join('、')+' 键'+(flip?'（翻环后）':'');
  const nAxial=types.filter(x=>x==='a').length;
  const penalty=subs.reduce((sum,s)=>sum+(slotAt(s.k,s.face,flip)==='a'?AVAL[s.kind]:0),0);
  document.getElementById('disub-na').innerHTML=nAxial+' <small>个</small>';
  document.getElementById('disub-energy').innerHTML=penalty.toFixed(1)+' <small>kJ/mol</small>';
  document.querySelectorAll('#g1-row [data-g]').forEach(b=>b.classList.toggle('on',b.dataset.g===disub.g1));
  document.querySelectorAll('#g2-row [data-g]').forEach(b=>b.classList.toggle('on',b.dataset.g===disub.g2));
  let pills;
  if(disub.t<.05||disub.t>.95){
    if(nAxial===0)pills='<span class="strain-pill cool">ee · 双 e 键优势构象</span>';
    else if(nAxial===2)pills='<span class="strain-pill hot">aa · 双 a 键，翻环后得 ee</span>';
    else pills='<span class="strain-pill warm">ea/ae · 总有一个直立基团</span>';
    // 最大基团是否被锁定在 e 键
    let big=subs[0]; subs.forEach(s=>{if(AVAL[s.kind]>AVAL[big.kind])big=s;});
    if(AVAL[big.kind]>=9){
      const ty=slotAt(big.k,big.face,flip);
      pills+=ty==='e'
        ?`<span class="strain-pill cool">${GROUP_CN[big.kind]}锁定 e 键</span>`
        :`<span class="strain-pill hot">${GROUP_CN[big.kind]}在 a 键（极少存在）</span>`;
    }
  }else pills='<span class="strain-pill warm">翻环中：顺反不变，a ⇄ e 互换</span>';
  document.getElementById('disub-strains').innerHTML=pills;
}
function flipDisub(){
  stopTweens();const start=disub.t,end=disub.t>.5?0:1;
  tweenTo(900,t=>{disub.t=lerp(start,end,t);renderDisub();});
}
// 构建两排基团选择 chip(显示基团符号 + A 值)
function buildGroupRows(){
  ['g1-row','g2-row'].forEach((rid,row)=>{
    const wrap=document.getElementById(rid);
    wrap.innerHTML=DISUB_GROUPS.map(g=>
      `<button class="chip${(row?disub.g2:disub.g1)===g?' on':''}" data-g="${g}" title="${GROUP_CN[g]} · A=${AVAL[g]} kJ/mol">${GROUP_LABEL[g]}<small> ${AVAL[g]}</small></button>`
    ).join('');
    wrap.querySelectorAll('[data-g]').forEach(b=>b.addEventListener('click',()=>{
      if(row===0)disub.g1=b.dataset.g; else disub.g2=b.dataset.g;
      disub.t=0; renderDisub();
    }));
  });
}
function initDisub(){
  buildGroupRows();
  bindToolbar('tb-disub',viewers.disub,act=>{if(act==='flip')flipDisub();});
  document.querySelectorAll('[data-pos12]').forEach(b=>b.addEventListener('click',()=>{
    disub.pos=b.dataset.pos12;disub.t=0;
    document.querySelectorAll('[data-pos12]').forEach(x=>x.classList.toggle('on',x===b));
    renderDisub();
  }));
  document.querySelectorAll('#disub-ct-seg button').forEach(b=>b.addEventListener('click',()=>{
    disub.ct=b.dataset.ct;disub.t=0;
    document.querySelectorAll('#disub-ct-seg button').forEach(x=>x.classList.toggle('on',x===b));
    renderDisub();
  }));
  document.getElementById('load-tbu').addEventListener('click',()=>{
    disub.g1='Me';disub.g2='tBu';disub.pos='1,4';disub.ct='cis';disub.t=0;
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
  document.body.dataset.page=name;   // 供 CSS 按页微调(如手机端 06 页更紧凑)
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
  // 移动端操作提示:让用户知道可以单指拖拽旋转模型
  document.querySelectorAll('.hint3d').forEach(h=>{
    h.textContent = isMobile ? '单指旋转 · 双指拖动平移/捏合缩放 · 双击复位' : '左键旋转 · 右键平移 · 滚轮缩放 · 双击复位';
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

/* ============================================================
   中英双语切换
============================================================ */
const I18N={
  nav:{intro:['原理','Basics'],ethane:['乙烷','Ethane'],butane:['丁烷','Butane'],
    rings:['环烷烃','Rings'],chair:['椅/船','Chair/Boat'],axeq:['a/e键','a/e Bonds'],disub:['二取代','Disub.']},
  navsub:{intro:[' · 核心',' · Core'],ethane:[' · 扭转张力',' · Torsion'],butane:[' · 空间位阻',' · Steric'],
    rings:[' · 角张力',' · Angle'],chair:[' · 翻环能垒',' · Ring Flip'],axeq:[' · 取代基',' · Substituent'],disub:[' · 顺反',' · cis/trans']},
  brand:'烷烃与环烷烃 · 空间结构互动课|Alkanes & Cycloalkanes · Interactive',
  toggle:{vdw:['◯ 填充:关','◯ Fill:Off'],pan:['✥ 平移','✥ Pan'],panOn:['✥ 平移中','✥ Panning'],
    angle:['∠ 键角','∠ Angles'],clash:['⚠ 冲突线','⚠ Clash Lines'],label:['🏷 标注','🏷 Labels'],
    rot:['⟳ 视角自转','⟳ Auto-Rotate'],flag:['⚑ 高亮旗杆氢','⚑ Flag H'],aecolor:['🎨 a/e 着色','🎨 a/e Color'],
    reset:['复位视角','Reset View'],flip:['🔄 翻环','🔄 Flip'],
    understand:['我已理解本节','I Understand This']},
  // 标题/标签短语:完整匹配(中文原文 → 英文)
  phrases:{
    '控制台':'Console','原理讲解':'Key Concepts','能量观点：优势构象如何决定？':'Energy: What Makes a Favored Conformation?',
    '规律总结（对应课件 P37–P38）':'Summary (Slides P37–P38)','NEWMAN 投影式（沿 C–C 键看）':'NEWMAN PROJECTION (along C–C)',
    'NEWMAN 投影式（沿 C2–C3 键看）':'NEWMAN PROJECTION (along C2–C3)','能量曲线（点击曲线可跳转）':'Energy Curve (click to seek)',
    '扭转角':'Torsion','二面角':'Dihedral','当前状态':'Current State','构象名称':'Conformation','扭转张力能':'Torsional E',
    '取代位置':'Position','顺 / 反':'cis / trans','快捷示例':'Quick Example','直立基团数':'Axial groups','位阻惩罚':'Strain penalty',
    'C1 取代基 ':'C1 group ','平均键角':'Avg angle','最小键角':'Min angle','环大小':'Ring size','稳定性':'Stability',
    '键角偏离 109.5°':'Deviation from 109.5°','二面角 φ':'Dihedral φ','甲基距离':'Me–Me dist','能量':'Energy','最近 H···H':'Nearest H···H'
  }
};
function t(pair){ return LANG==='en'?pair[1]:pair[0]; }
function applyLang(){
  const en=LANG==='en';
  // 导航
  document.querySelectorAll('.nav-item[data-page]').forEach(b=>{
    const p=b.dataset.page; if(!I18N.nav[p])return;
    const sub=b.querySelector('.sub-t');
    b.childNodes.forEach(n=>{ if(n.nodeType===3 && n.textContent.trim() && !n._cn) n._cn=n.textContent; });
    const textNode=Array.from(b.childNodes).find(n=>n.nodeType===3 && n.textContent.trim());
    if(textNode) textNode.textContent=en?I18N.nav[p][1]:(textNode._cn||I18N.nav[p][0]);
    if(sub){ if(!sub._cn)sub._cn=sub.textContent; sub.textContent=en?I18N.navsub[p][1]:sub._cn; sub.style.display=''; }
  });
  // brand 副标题
  const bp=document.querySelector('.brand p');
  if(bp){ if(!bp._cn)bp._cn=bp.textContent; bp.textContent=en?I18N.brand.split('|')[1]:bp._cn; }
  // 工具栏按钮(按 data-toggle / data-act)
  const setBtn=(btn,cn,enTxt)=>{ if(!btn._cn)btn._cn=btn.textContent; btn.textContent=en?enTxt:btn._cn; };
  const TB2VW={'tb-ethane':'ethane','tb-butane':'butane','tb-rings':'rings','tb-chair':'chair','tb-axeq':'axeq','tb-disub':'disub'};
  document.querySelectorAll('[data-toggle]').forEach(b=>{
    const k=b.dataset.toggle;
    if(k==='vdw'){
      const tb=b.closest('.toolbar');
      const vw=tb&&viewers[TB2VW[tb.id]];
      if(vw) syncVdwBtn(b,vw);
      return;
    }
    if(k==='pan'){ setBtn(b,I18N.toggle.pan[0], b.classList.contains('on')?I18N.toggle.panOn[1]:I18N.toggle.pan[1]); return; }
    if(I18N.toggle[k]) setBtn(b,I18N.toggle[k][0],I18N.toggle[k][1]);
  });
  document.querySelectorAll('[data-act="resetView"]').forEach(b=>setBtn(b,I18N.toggle.reset[0],I18N.toggle.reset[1]));
  document.querySelectorAll('[data-act="flip"]').forEach(b=>setBtn(b,I18N.toggle.flip[0],I18N.toggle.flip[1]));
  document.querySelectorAll('.ubtn').forEach(b=>setBtn(b,I18N.toggle.understand[0],I18N.toggle.understand[1]));
  // 短语扫描(card-h / h4 / 表头 / stat 标签)
  document.querySelectorAll('.card-h, .ctl-block h4, .stat .k, table th, .state-en, .k, .range-row span').forEach(el=>{
    const raw=el.textContent.trim();
    if(I18N.phrases[raw]!==undefined){
      if(!el._cn)el._cn=raw;
      // card-h 含 dot span,保留元素子节点,只替换文字
      const dot=el.querySelector('.dot');
      el.textContent=en?I18N.phrases[raw]:el._cn;
      if(dot)el.insertBefore(dot,el.firstChild);
    }
  });
  const lb=document.getElementById('lang-btn'); if(lb)lb.textContent=en?'EN':'中';
}
(function initLang(){
  LANG=localStorage.getItem('conformation-lab-lang')||'cn';
  const lb=document.getElementById('lang-btn');
  if(lb)lb.addEventListener('click',()=>{
    LANG=LANG==='en'?'cn':'en';
    localStorage.setItem('conformation-lab-lang',LANG);
    // 重建动态内容后再套用界面语言
    try{renderEthane();renderButane();renderRings();renderChair();renderAxeq();renderDisub();}catch(e){}
    applyLang();
  });
  if(LANG==='en')applyLang();   // 记住英文偏好:首屏即英文
})();

/* ============================================================
   日间/夜间模式切换
============================================================ */
const THEME_KEY='conformation-lab-theme';
function applyTheme(light){
  document.body.classList.toggle('light',light);
  const btn=document.getElementById('theme-btn');
  if(btn) btn.textContent=light?'☀️':'🌙';
  // 同步更新所有 viewer 的场景背景
  Object.values(viewers).forEach(v=>v._updateBg && v._updateBg());
}
(function initTheme(){
  // 默认跟随系统,已保存的优先
  const saved=localStorage.getItem(THEME_KEY);
  const light = saved!==null ? saved==='light' : matchMedia('(prefers-color-scheme: light)').matches;
  applyTheme(light);
  const btn=document.getElementById('theme-btn');
  if(btn) btn.addEventListener('click',()=>{
    const now=!document.body.classList.contains('light');
    localStorage.setItem(THEME_KEY,now?'light':'dark');
    applyTheme(now);
  });
})();

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
