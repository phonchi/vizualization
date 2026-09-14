/* Earthquake science exhibits. Fixed datasets; interaction changes the view,
   never fits or recalibrates a model. No network, dependencies, or inline scripts. */
(function (global) {
'use strict';
// Notebook outputs may embed the shared runtime beside the global site asset.
// Reuse the first browser controller; Node rendering/export stays independent.
if(typeof document!=='undefined' && global.EarthquakeExhibits){
 global.EarthquakeExhibits.mount(document);
 return;
}
const esc=v=>String(v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const clamp=(x,a=0,b=1)=>Math.max(a,Math.min(b,x));
const fmt=(x,n=2)=>Number(x).toFixed(n);
// Native range step is 0.1 percent: preserve it when reflecting state to DOM.
const progressValue=p=>Math.round(clamp(p)*1000)/10;
const txt=(x,y,t,c='ex-label',anchor='start')=>`<text x="${x}" y="${y}" class="${c}" text-anchor="${anchor}">${esc(t)}</text>`;
const line=(x1,y1,x2,y2,c='ex-axis')=>`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" class="${c}"/>`;
const rect=(x,y,w,h,c='ex-tile',extra='')=>`<rect x="${x}" y="${y}" width="${Math.max(0,w)}" height="${Math.max(0,h)}" rx="3" class="${c}" ${extra}/>`;
const circle=(x,y,r,c='ex-event',extra='')=>`<circle cx="${x}" cy="${y}" r="${r}" class="${c}" ${extra}/>`;
const path=(pts,c='ex-curve')=>`<path d="${pts.map((p,i)=>(i?'L':'M')+p.map(v=>Number(v).toFixed(2)).join(' ')).join(' ')}" class="${c}"/>`;
const curve=(f,x0,x1,y0,scale,n=100,c='ex-curve')=>path(Array.from({length:n+1},(_,i)=>{const t=i/n;return[x0+t*(x1-x0),y0-scale*f(t)];}),c);
const flow=(a,b,c='ex-link')=>`<path d="M${a[0]},${a[1]} C${a[0]},${(a[1]+b[1])/2} ${b[0]},${(a[1]+b[1])/2} ${b[0]},${b[1]}" class="${c}"/>`;
const grid=(w,h,x=44,y=24)=>[0,.25,.5,.75,1].map(t=>line(x,y+t*h,x+w,y+t*h,'ex-grid')).join('');
const labelNumber=(x,y,num,label)=>txt(x,y,num,'ex-number')+txt(x,y+25,label,'ex-small');
const HAWKES_EVENTS=[2,5,5.8];
function hawkesAt(t){
 const history=HAWKES_EVENTS.filter(e=>e<t);
 const components=history.map(e=>1.5*Math.exp(-(t-e)/1.3));
 return{background:.3,history,components,total:.3+components.reduce((a,b)=>a+b,0)};
}
const range=Array.from;
const kinds={
 pipeline:{play:true,label:'敘事步驟',start:.65},map:{play:false,label:'歷史截止年',start:1,layers:['S','R','events']},timeline:{play:true,label:'發報窗',start:.46},
 magnitude:{play:true,label:'原始規模',start:.49},integration:{play:true,label:'累積進度',start:.65},distribution:{play:false,label:'觀測值',start:.333,select:['下尾','上尾','雙尾']},
 kernel:{play:true,label:'剖面位置',start:.5,layers:['contours','sources']},hawkes:{play:true,label:'时间（天）',start:.75},branching:{play:true,label:'世代',start:.7},
 product:{play:true,label:'三核的組合',start:.7},matrix:{play:true,label:'模擬目錄',start:.55,select:['N','S','M','cL','L']},mix:{play:true,label:'固定混合的形成',start:.6},
 family:{play:false,label:'模型焦點',start:.3,select:['SUP','PPE','ETAS','EEPAS','混合']},protocol:{play:true,label:'評估時間位置',start:.6,select:['共同時鐘','延遲驗證','部分前瞻']}
};
function render(kind,p,s={},data={},w=760){
 p=clamp(p);w=Math.max(280,w);const mobile=w<600;let h=360,body='',caption='',value='';const left=44,right=w-28,cw=right-left;
 if(kind==='hawkes'){
  h=mobile?410:360;const end=8,t=p*end,baseline=h-105,top=40,scale=(baseline-top)/4.3,events=HAWKES_EVENTS;
  body+=grid(cw,baseline-top,left,top)+line(left,baseline,right,baseline)+txt(left,22,'條件率 λ*(t)','ex-small');
  const X=x=>left+x/end*cw;
  body+=rect(X(t),top,right-X(t),baseline-top,'ex-future');
  body+=line(left,baseline-.3*scale,right,baseline-.3*scale,'ex-dashed');
  for(let i=0;i<events.length;i++)if(events[i]<t){const e=events[i];const pts=[];for(let x=e;x<=t;x+=.025)pts.push([X(x),baseline-1.5*Math.exp(-(x-e)/1.3)*scale]);body+=path(pts,'ex-component ex-component-'+i);}
  const points=[];for(let x=0;x<=t;x+=.015){points.push([X(x),baseline-hawkesAt(x).total*scale]);}body+=path(points,'ex-curve ex-strong');
  body+=line(X(t),top,X(t),baseline+40,'ex-cursor');
  events.forEach((e,i)=>{if(e<t)body+=circle(X(e),baseline+36,6,'ex-event')+txt(X(e),baseline+64,String(e),'ex-small','middle');});
  body+=txt(left,h-8,'已知歷史 Hₜ','ex-small')+txt(right,h-8,`t = ${fmt(t)} 天`,'ex-small','end');
  const rate=hawkesAt(t).total;value=fmt(rate)+' 顆／天';
  caption=t<2?'尚無事件：只有背景率 0.3。':`${events.filter(e=>e<t).length} 顆事件已進入歷史。細線各是一顆事件的衰減貢獻；虛線為背景 0.3，粗線將兩者相加。歷史只含嚴格早於 t 的事件。`;
 } else if(kind==='map'){
  h=mobile?480:490;const pts=data.polygon||[],cells=data.cells||[],events=data.events||[];const xs=pts.map(p=>p[0]),ys=pts.map(p=>p[1]);const minx=Math.min(...xs),maxx=Math.max(...xs),miny=Math.min(...ys),maxy=Math.max(...ys);const scale=Math.min((w-56)/(maxx-minx),(h-70)/(maxy-miny));const ox=(w-(maxx-minx)*scale)/2,oy=26;const P=([x,y])=>[ox+(x-minx)*scale,oy+(maxy-y)*scale];
  const year=Math.round(1960+p*61);body+=txt(24,24,'義大利 · 投影公里','ex-small');
  if(s.S!==false)body+=path([...pts,pts[0]].map(P),'ex-region-s');
  if(s.R!==false)cells.forEach((c,i)=>{const a=P([c[0],c[3]]),b=P([c[1],c[2]]);body+=rect(a[0],a[1],b[0]-a[0],b[1]-a[1],'ex-region-r',`data-cell="${i}"`);});
  const visible=events.filter(e=>e[2]<=year);if(s.events!==false)visible.forEach(e=>{const xy=P(e);body+=circle(...xy,2.4+(e[3]-5)*2,e[5]>=0&&e[2]>=2012?'ex-event':'ex-source');});
  const bar=100*scale;body+=line(24,h-26,24+bar,h-26,'ex-scale')+txt(24,h-38,'100 km','ex-small');
  value=year+' 年';caption=`截至 ${year}：${visible.length} 顆 M≥5 歷史事件。藍格是實際 177 格，外框是 CPTI15 收集區；框外未顯示海岸線。`;
 } else if(kind==='timeline'){
  h=mobile?430:350;const win=Math.min(39,Math.floor(p*40)),datawin=(data.windows||[])[win];const X=i=>left+i/40*cw;
  body+=txt(left,26,'2012','ex-small')+txt(right,26,'2021','ex-small','end');
  for(let i=0;i<40;i++)body+=rect(X(i),44,cw/40-1,38,i===win?'ex-active':'ex-window');
  (data.targets||[]).forEach(e=>{const q=(Date.parse(e[4].replace(' ','T'))-Date.UTC(2012,0,1))/86400000/91.31;if(q>=0&&q<40)body+=circle(X(q),105+(e[3]-5)*-8,4,'ex-event');});
  const y=mobile?172:155,cut=left+cw*.33,issue=left+cw*.58;
  body+=txt(left,y-20,`第 ${win+1} 窗：放大資料截止`,'ex-label');body+=rect(left,y,issue-left,52,'ex-history')+rect(cut,y,issue-cut,52,'ex-delay')+rect(issue,y,right-issue,52,'ex-active');
  body+=line(cut,y-8,cut,y+68,'ex-dashed')+line(issue,y-8,issue,y+68,'ex-cursor');
  body+=txt(cut,y+94,'提前 50 天','ex-small','middle')+txt(issue,y+122,'發報','ex-label','middle')+txt(right,y+94,'91.31 天後','ex-small','end');
  value='第 '+(win+1)+'／40 窗';caption=`${datawin?datawin.start.slice(0,10):''} 發報。PPE／EEPAS 不用灰褐色延遲帶內的來源；ETAS 可使用發報前历史。红點是事後目標，不能提前看見。`.replace('历史','歷史').replace('红','紅');
 } else if(kind==='magnitude'){
  h=mobile?370:310;const raw=+(2.2+p*3.6).toFixed(2),binned=Math.floor(raw*10+.5+1e-8)/10;const origin=Math.floor(raw*10)/10-.2;const y=85,bw=cw/7;
  body+=txt(left,24,'原始規模 → 十分位取整','ex-label');for(let i=0;i<7;i++){const c=+(origin+i*.1).toFixed(1),active=Math.abs(c-binned)<.001;body+=rect(left+i*bw,y,bw-2,86,active?'ex-active':'ex-window')+txt(left+(i+.5)*bw,y+54,c.toFixed(1),'ex-small','middle');}
  const pos=left+(raw-(origin-.05))/.7*cw;body+=line(pos,y-22,pos,y-2,'ex-cursor')+circle(pos,y-24,5,'ex-event');
  body+=labelNumber(left,232,raw.toFixed(2),'原始 Mw')+labelNumber(mobile?w*.55:w*.42,232,binned.toFixed(1),'取整規模');
  value=binned>=5?'進入目標集合':binned>=2.5?'可作輸入事件':'低於輸入門檻';caption=`${raw.toFixed(2)} 取整為 ${binned.toFixed(1)}。原始 2.45／4.95 對應名目 2.5／5.0；圖中只移動觀測值，不更改實驗門檻。`;
 } else if(kind==='integration'){
  h=mobile?560:360;const size=Math.min(mobile?w-70:w*.44,290),x=mobile?(w-size)/2:30,y=26,n=6,step=size/n,filled=Math.round(p*n*n);let total=0,all=0;
  for(let j=0;j<n;j++)for(let i=0;i<n;i++){const v=.006+.04*Math.exp(-((i-2)**2+(j-3)**2)/3);all+=v;const idx=j*n+i;if(idx<filled)total+=v;body+=rect(x+i*step,y+j*step,step-3,step-3,idx<filled?'ex-active':'ex-window',`style="fill-opacity:${idx<filled?.3+v*13:.2}"`);}
  const xx=mobile?30:w*.56,yy=mobile?size+82:96;body+=labelNumber(xx,yy,fmt(total,3),'目前累積的期望數 Λ');
  const base=yy+143,bwid=(mobile?w-70:w*.4)/4;for(let k=0;k<4;k++){let pr=Math.exp(-total)*Math.pow(total,k);for(let i=2;i<=k;i++)pr/=i;body+=rect(xx+k*bwid,base-pr*85,bwid*.6,pr*85,'ex-bar')+txt(xx+k*bwid+bwid*.3,base+25,k+' 顆','ex-small','middle');}
  value=filled+'／36 小格';caption=`每塊面積乘上其率密度，再累積。完整固定示例 Λ=${fmt(all,3)}；右側是採 Poisson 假設後的 0–3 顆機率。`;
 } else if(kind==='distribution'){
  h=mobile?430:360;const counts=[1,2,4,6,4,2,1],obs=Math.round(p*6),upper=s.variant==='上尾',two=s.variant==='雙尾',bw=cw/7,base=180,max=6;
  body+=txt(left,24,'模擬統計量的分布','ex-label');let accum=0,cdf=[];counts.forEach((c,i)=>{const isTail=two?(i<=Math.min(obs,6-obs)||i>=Math.max(obs,6-obs)):upper?i>=obs:i<=obs;const xx=left+i*bw;body+=rect(xx+3,base-c/max*110,bw-6,c/max*110,isTail?'ex-active':'ex-window')+txt(xx+bw/2,base+24,i,'ex-small','middle');accum+=c;cdf.push([xx+bw/2,330-accum/20*90]);});
  body+=line(left,base,right,base)+txt(left,240,'累積機率 CDF','ex-small');const stair=[[left,330]];cdf.forEach(pt=>{stair.push([pt[0],stair[stair.length-1][1]],pt);});stair.push([right,cdf[cdf.length-1][1]]);body+=path(stair,'ex-curve')+circle(...cdf[obs],6,'ex-focus-dot');
  const q=counts.slice(0,obs+1).reduce((a,b)=>a+b,0)/20,u=counts.slice(obs).reduce((a,b)=>a+b,0)/20;value=`${s.variant||'下尾'} = ${fmt(two?Math.min(1,2*Math.min(q,u)):upper?u:q)}`;caption=`觀測值 ${obs}，q=${fmt(q)}；同值比例 ${fmt(counts[obs]/20)}。上尾=${fmt(u)}，包含等號。雙尾採 2×較小尾端、上限為 1。`;
 } else if(kind==='kernel'){
  h=mobile?560:370;const size=Math.min(mobile?w-48:w*.46,310),ox=24,oy=30,cut=p;const centers=[[.30,.38],[.70,.66]];
  body+=txt(ox,20,'固定兩個來源的等值圈','ex-small');if(s.contours!==false)centers.forEach((c,i)=>[.10,.18,.26].forEach(r=>body+=`<ellipse cx="${ox+c[0]*size}" cy="${oy+c[1]*size}" rx="${r*size}" ry="${r*size}" class="ex-contour ex-component-${i}"/>`));
  if(s.sources!==false)centers.forEach(c=>body+=circle(ox+c[0]*size,oy+c[1]*size,6,'ex-source'));
  body+=line(ox,oy+cut*size,ox+size,oy+cut*size,'ex-cursor');
  const x0=mobile?44:w*.56,x1=right,base=mobile?size+205:285,scale=100;body+=txt(x0,base-122,'同一截線上的貢獻','ex-small')+line(x0,base,x1,base);
  const fn=(u,i)=>Math.exp(-((u-centers[i][0])**2+(cut-centers[i][1])**2)/.035);
  centers.forEach((_,i)=>body+=curve(u=>fn(u,i),x0,x1,base,scale,80,'ex-component ex-component-'+i));body+=curve(u=>fn(u,0)+fn(u,1),x0,x1,base,scale,80);
  value='剖面 '+Math.round(p*100)+'%';caption='移動的是閱讀截線。兩個核及其參數固定；粗線是逐點總和，細線保留各來源。此 Gaussian 示意不代替義大利 PPE 核。';
 } else if(kind==='branching'){
  h=mobile?490:430;const gen=Math.min(3,Math.floor(p*4)),nodes=[{x:.5,y:45,g:0,parent:-1},{x:.25,y:145,g:1,parent:0},{x:.75,y:145,g:1,parent:0},{x:.12,y:245,g:2,parent:1},{x:.38,y:245,g:2,parent:1},{x:.75,y:245,g:2,parent:2},{x:.38,y:345,g:3,parent:4}];
  nodes.forEach((n,i)=>{if(n.g<=gen&&n.parent>=0){const a=nodes[n.parent];body+=flow([left+a.x*cw,a.y+14],[left+n.x*cw,n.y-14]);}});
  nodes.forEach((n,i)=>{if(n.g<=gen)body+=circle(left+n.x*cw,n.y,i?11:17,i?'ex-event':'ex-focus-dot')+txt(left+n.x*cw,n.y+36,i?'事件 '+i:'背景事件','ex-small','middle');});
  value='第 '+gen+' 代';caption=`目前顯示 ${nodes.filter(n=>n.g<=gen).length} 個固定示意事件。連線是生成時指定的親代關係；這棵樹的後代數不是分支比估計。`;
 } else if(kind==='product'){
  h=mobile?650:370;const phase=Math.min(3,Math.floor(p*4)),panels=mobile?[[24,24,w-48],[24,190,w-48],[24,356,w-48]]:[[24,26,w/3-32],[w/3+8,26,w/3-32],[2*w/3-8,26,w/3-32]];
  const titles=['規模核 g(m)','時間核 f(τ)','空間核 h(x,y)'];panels.forEach(([x,y,ww],i)=>{body+=txt(x,y+20,titles[i],'ex-label');const base=y+126;body+=line(x,base,x+ww,base);if(i<2)body+=curve(u=>i===0?Math.exp(-.5*((u-.52)/.14)**2):(u===0?0:Math.exp(-.5*((Math.log(u)+1.2)/.45)**2)/u/4),x,x+ww,base,78,90,i<=phase?'ex-curve':'ex-dashed');else [20,36,52].forEach(r=>body+=`<circle cx="${x+ww/2}" cy="${y+82}" r="${r}" class="${i<=phase?'ex-contour':'ex-dashed'}"/>`);if(i===phase)body+=rect(x,y+145,ww,4,'ex-active');});
  const yy=mobile?565:260;body+=txt(w/2,yy,phase<3?['先分配未來規模','再分配等待時間','再分配空間位置'][phase]:'單顆貢獻 ∝ g × f × h','ex-label','middle');
  if(phase===3)for(let j=0;j<3;j++)for(let i=0;i<9;i++)body+=rect(w/2-108+i*24,yy+22+j*16,20,12,'ex-active',`style="fill-opacity:${.15+.7*Math.exp(-((i-4)**2+(j-1)**2)/4)}"`);
  value=['規模','時間','空間','乘積'][phase];caption='一個固定輸入事件，提供三個方向的密度。三核相乘後仍須乘權重、積分到格箱，再加其他事件與背景；不是三份獨立預報。';
 } else if(kind==='matrix'){
  h=mobile?430:360;const tests=['N','S','M','cL','L'],chosen=tests.includes(s.variant)?s.variant:'N',cols=['總數','位置','規模','時間'],rules={N:['抽樣','—','—','合計'],S:['固定','抽樣','合計','合計'],M:['固定','合計','抽樣','合計'],cL:['固定','抽樣','抽樣','抽樣'],L:['抽樣','抽樣','抽樣','抽樣']};const labelw=38,cellw=(cw-labelw)/4;
  cols.forEach((c,i)=>body+=txt(left+labelw+(i+.5)*cellw,30,c,'ex-small','middle'));tests.forEach((t,j)=>{const y=52+j*43;body+=txt(left,y+26,t,'ex-label');rules[t].forEach((r,i)=>{body+=rect(left+labelw+i*cellw,y,cellw-3,37,t===chosen?'ex-active-soft':'ex-window')+txt(left+labelw+(i+.5)*cellw,y+25,r,'ex-small','middle');});});
  const count=chosen==='N'||chosen==='L'?[6,4,8,5][Math.min(3,Math.floor(p*4))]:6,y=mobile?350:310;for(let i=0;i<count;i++){const bin=(i*3+Math.floor(p*12))%8;body+=circle(left+(bin+.5)/8*cw,y+(i%2)*20,5,'ex-event');}
  value=chosen+'-test';caption=`${chosen}：${rules[chosen].map((r,i)=>cols[i]+r).join('、')}。下方點列用固定示意目錄展示條件；沒有用它重算義大利分數。`;
 } else if(kind==='mix'){
  h=mobile?480:360;const A=[2,6,1,3],B=[6,2,5,3],weight=.5,stage=p,base=mobile?325:260,bw=cw/4,scale=24;
  body+=txt(left,28,'固定权重：A ½ ＋ B ½'.replace('权','權'),'ex-label');for(let i=0;i<4;i++){const av=A[i]*weight,bv=B[i]*weight;body+=rect(left+i*bw+12,base-av*scale,bw-24,av*scale,'ex-mix-a');body+=rect(left+i*bw+12,base-av*scale-bv*scale*stage,bw-24,bv*scale*stage,'ex-mix-b');body+=txt(left+(i+.5)*bw,base+28,'格 '+(i+1),'ex-small','middle')+txt(left+(i+.5)*bw,base-(av+bv*stage)*scale-12,fmt(av+bv*stage,1),'ex-label','middle');}
  body+=line(left,base,right,base);value='Λ = '+fmt(6+8*p,1);caption=p===1?'固定一半混合：總期望數 ½×12＋½×16＝14。不是把兩份對數概似相加。':'播放正把 B 的一半貢獻疊到 A 的一半上。終點權重固定 ½；過程中的畫面只是加法演示。';
 } else if(kind==='family'){
  h=mobile?490:370;const names=['SUP','PPE','ETAS','EEPAS','混合'],focus=names[Math.min(4,Math.floor(p*5))],pos=mobile?[[.2,65],[.75,65],[.2,220],[.75,220],[.5,390]]:[[.12,90],[.38,65],[.64,90],[.88,65],[.5,285]];const links=[[0,1],[0,2],[1,3],[2,4],[3,4]];
  links.forEach(([i,j])=>{const a=pos[i],b=pos[j];body+=flow([a[0]*w,a[1]+22],[b[0]*w,b[1]-22],names[i]===focus||names[j]===focus?'ex-link ex-emphasis':'ex-link');});
  pos.forEach((p,i)=>{body+=circle(p[0]*w,p[1],names[i]===focus?34:25,names[i]===focus?'ex-focus-ring':'ex-model-ring')+txt(p[0]*w,p[1]+6,names[i],'ex-label','middle');});
  value=focus;caption={SUP:'空間均勻、固定學習率。它提供一個所有模型共用的簡單比較起點。',PPE:'依來源目錄缓慢更新空間權重，本站採50天資料延遲。'.replace('缓','緩'),ETAS:'將背景與事件後的觸發貢獻相加；本站預報使用已知歷史第一代近似。',EEPAS:'沿規模、時間與空間的前兆尺度分配貢獻；PPE提供背景形狀。',混合:'結合ETAS與EEPAS的不同資訊；成分和權重須在評估前指定。'}[focus];
 } else if(kind==='protocol'){
  h=mobile?530:380;
  const partial=s.variant==='部分前瞻',delayed=s.variant==='延遲驗證';
  const names=['回溯','擬前瞻',partial?'部分前瞻':'前瞻'],top=50,step=mobile?135:92,start=mobile?70:120,end=w-28,xx=start+p*(end-start);
  body+=txt(start,25,'過去 → 未來','ex-small');
  names.forEach((n,i)=>{
   const y=top+i*step;
   body+=partial&&i===2&&mobile?txt(12,y+18,'部分','ex-label')+txt(12,y+40,'前瞻','ex-label'):txt(12,y+25,n,'ex-label');
   body+=line(start,y+25,end,y+25);
   const marker=start+(end-start)*(i===0?.9:i===1||partial?.42:.16);
   if(i===2&&partial){
    body+=rect(start,y+9,marker-start,32,'ex-history')+rect(marker,y+9,end-marker,32,'ex-active-soft');
    body+=txt((start+marker)/2,y+30,'已知','ex-small','middle')+txt((marker+end)/2,y+30,'未知','ex-small','middle');
   }else{
    body+=rect(start,y+9,(i===0?end:xx)-start,32,i===0?'ex-history':'ex-active-soft');
    if(i===1)body+=rect(xx,y+9,end-xx,32,'ex-future');
    if(i===2&&delayed)body+=rect(Math.max(start,xx-(end-start)*.2),y+9,Math.min(xx-start,(end-start)*.2),32,'ex-delay');
   }
   body+=line(marker,y-2,marker,y+45,'ex-cursor');
   const note=i===0?'評估既有資料':i===1?'依事件截止重演':partial?'預報涵蓋已知與未知':delayed?'模型與建模者隔離':'先鎖定流程';
   body+=txt(marker,y+68,note,'ex-small','middle');
  });
  body+=line(xx,top-10,xx,top+step*2+50,'ex-current-time');value='評估位置 '+Math.round(p*100)+'%';
  caption=partial?'部分前瞻同時評估已知與未知時段。例如主震後幾天才發報，評分卻從主震時刻開始；已知幾天與真正未來的部分須分開說明。詳見本章的檢驗方式對照。':delayed?'受控延遲前瞻：事先鎖定流程與測試規格，測試資料對模型及建模者保持受控隔離，再由受控程序稍後執行。一般晚點評分不等於延遲前瞻。':'回溯是評估既有資料，可以是樣本內或樣本外；擬前瞻是其中依事件截止重演的子類。本站用修訂 HORUS 回放事件時間，並未逐版還原歷史當時可得的目錄。';
 } else if(kind==='pipeline'){
  h=mobile?570:380;const stages=['讀目錄','訂規格','模型','發預報','檢驗'],active=Math.min(4,Math.floor(p*5));const positions=mobile?[[w/2,35],[w/2,125],[w/2,255],[w/2,390],[w/2,500]]:[[60,165],[w*.27,165],[w*.5,165],[w*.73,165],[w-55,165]];
  for(let i of [0,3])body+=flow(positions[i],positions[i+1],i<active?'ex-link ex-emphasis':'ex-link');
  const model=positions[2],offset=mobile?68:62;['SUP','PPE','ETAS','EEPAS'].forEach((n,i)=>{const x=model[0]+(i%2?1:-1)*offset,y=model[1]+(i<2?-44:44);body+=flow(positions[1],[x,y-25],active>=2?'ex-link ex-emphasis':'ex-link')+flow([x,y+25],positions[3],active>=3?'ex-link ex-emphasis':'ex-link');body+=circle(x,y,25,active>=2?'ex-model-ring':'ex-window')+txt(x,y+5,n,'ex-small','middle');});
  positions.forEach((q,i)=>{if(i!==2)body+=circle(q[0],q[1],active===i?30:24,active>=i?'ex-focus-ring':'ex-model-ring')+txt(q[0],q[1]+5,stages[i],'ex-small','middle');});
  value=stages[active];caption=['先讀資料來源、規模尺度與收錄能力。','先定區域、門檻、發報窗和評分方法。','四個模型回答同一份規格；形狀與歷史用法各有不同。','每格、每箱寫下期望數，保留發報時刻。','用之後的目標事件核對數量、分配與相對資訊。'][active];
 }
 return{body,height:h,caption,value};
}
function svgHTML(kind,p,state,data,w=760,uid='exhibit'){
 const r=render(kind,p,state,data,w);return{...r,svg:`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${r.height}" role="img" aria-labelledby="${uid}-title ${uid}-desc"><title id="${uid}-title">${esc(state.title||kind)}</title><desc id="${uid}-desc">${esc(r.caption)}</desc>${r.body}</svg>`};
}
const instances=new Set();let serial=0;
function mount(scope){if(typeof document==='undefined')return;scope=scope||document;const roots=scope.matches&&scope.matches('.quake-exhibit')?[scope]:[...scope.querySelectorAll('.quake-exhibit')];roots.forEach(root=>{
 if(root.dataset.mounted)return;root.dataset.mounted='true';const config=JSON.parse(root.dataset.config||'{}'),kind=root.dataset.kind,def=kinds[kind];if(!def)return;
 const state={title:config.title,variant:def.select?def.select[0]:undefined},stage=root.querySelector('.ex-stage'),slider=root.querySelector('[data-control="progress"]'),readout=root.querySelector('.ex-readout'),caption=root.querySelector('.ex-caption'),play=root.querySelector('[data-action="play"]');
 let p=Number(slider.value)/100,playing=false,raf=0,last=0,visible=true;const uid=root.dataset.exhibit+'-'+(++serial),reduced=matchMedia('(prefers-reduced-motion: reduce)');
 function draw(announce=false){const width=Math.max(280,Math.round(stage.clientWidth||760)),r=svgHTML(kind,p,state,config.data||{},width,uid);stage.innerHTML=r.svg;readout.textContent=r.value;caption.textContent=r.caption;slider.value=String(progressValue(p));slider.setAttribute('aria-valuetext',r.value);if(kind==='family'){const q=root.querySelector('[data-control=compare]');if(q)q.value=r.value;}root.dataset.progress=fmt(p,3);if(announce)root.querySelector('.ex-status').textContent=r.value;}
 function stop(){playing=false;cancelAnimationFrame(raf);raf=0;last=0;if(play){play.textContent=reduced.matches?'下一個畫面':'播放導覽';play.setAttribute('aria-pressed','false');}root.dataset.playing='false';}
 function tick(time){if(!playing)return;if(!visible||document.hidden){stop();return;}if(last)p=clamp(p+(time-last)/14000);last=time;draw();if(p>=1){stop();return;}raf=requestAnimationFrame(tick);}
 if(play)play.addEventListener('click',()=>{if(playing){stop();return;}if(reduced.matches){p=p>=1?0:clamp(p+.1);draw(true);return;}if(p>=1)p=0;playing=true;play.textContent='暫停';play.setAttribute('aria-pressed','true');root.dataset.playing='true';raf=requestAnimationFrame(tick);});
 slider.addEventListener('input',()=>{stop();p=Number(slider.value)/100;draw(true);});root.querySelector('[data-action="reset"]').addEventListener('click',()=>{stop();p=0;draw(true);});
 root.querySelectorAll('[data-layer]').forEach(el=>{state[el.dataset.layer]=el.checked;el.addEventListener('change',()=>{state[el.dataset.layer]=el.checked;draw(true);});});
 root.querySelectorAll('[data-action=boundary]').forEach(el=>el.addEventListener('click',()=>{stop();p=(Number(el.dataset.value)-2.2)/3.6;draw(true);}));
 const select=root.querySelector('[data-control="compare"]');if(select)select.addEventListener('change',()=>{state.variant=select.value;if(kind==='family')p=def.select.indexOf(select.value)/4;draw(true);});
 const resize=new ResizeObserver(()=>draw());resize.observe(stage);const inter=new IntersectionObserver(entries=>{visible=entries[0].isIntersecting;if(!visible)stop();},{threshold:.05});inter.observe(root);
 const mode=()=>{stop();root.classList.toggle('ex-reduced',reduced.matches);if(reduced.matches)p=1;if(play)play.textContent=reduced.matches?'下一個畫面':'播放導覽';draw();};reduced.addEventListener('change',mode);mode();instances.add({root,stop,resize,inter});
 });}
function pauseAll(){instances.forEach(x=>x.stop());}
const api={mount,pauseAll,render,svgHTML,kinds,progressValue,hawkesAt};global.EarthquakeExhibits=api;
if(typeof module!=='undefined'&&module.exports)module.exports=api;
if(typeof document!=='undefined'){
 const boot=()=>{mount(document);new MutationObserver(records=>{for(const r of records)for(const n of r.addedNodes)if(n.nodeType===1&&!n.closest('.ex-stage'))mount(n);instances.forEach(x=>{if(!x.root.isConnected){x.stop();x.resize.disconnect();x.inter.disconnect();delete x.root.dataset.mounted;instances.delete(x);}});}).observe(document.body,{childList:true,subtree:true});};
 if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',boot);else boot();document.addEventListener('visibilitychange',()=>{if(document.hidden)pauseAll();});
}
})(typeof window!=='undefined'?window:globalThis);
