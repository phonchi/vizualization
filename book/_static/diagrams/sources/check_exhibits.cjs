const fs=require('fs'),path=require('path'),crypto=require('crypto');
const api=require('../../exhibits.js'),meta=JSON.parse(fs.readFileSync(path.join(__dirname,'exhibits.json'),'utf8'));const checks=[];
for(const m of meta){const def=api.kinds[m.kind];let states=0;for(const width of [320,390,760,1024])for(const p of [0,.2,.5,.8,1])for(const variant of def.select||[undefined]){
 const r=api.svgHTML(m.kind,p,{title:m.title,variant},m.data||{},width,'test-'+m.slug);
 if(/NaN|Infinity|undefined/.test(r.svg))throw Error(m.slug+' nonfinite output');
 if(!r.caption||!r.value||!r.body)throw Error(m.slug+' missing narrative');
 if(/(?:width|height)="-/.test(r.svg))throw Error(m.slug+' negative extent');
 if(!r.svg.includes('<title')||!r.svg.includes('<desc'))throw Error(m.slug+' absent accessible text');
 states++;
 }checks.push({slug:m.slug,states,passed:true});}
const map=meta.find(m=>m.kind==='map').data;if(map.cells.length!==177||map.events.length!==128||map.targets.length!==25||map.windows.length!==40)throw Error('Italy data cardinalities changed');
const h=api.render('hawkes',.125,{}, {},760);if(h.value!=='0.30 顆／天')throw Error('Hawkes pre-event background incorrect');
const mix=api.render('mix',1,{}, {},760);if(mix.value!=='Λ = 14.0')throw Error('Fixed mixture incorrect');
const d=api.render('distribution',1/3,{variant:'上尾'}, {},760);if(d.value!=='上尾 = 0.85')throw Error('Discrete upper tail incorrect');
const a=api.render('magnitude',(4.95-2.2)/3.6,{}, {},760);if(a.value!=='進入目標集合')throw Error('Magnitude boundary incorrect');
const b=api.render('magnitude',(2.45-2.2)/3.6,{}, {},760);if(b.value!=='可作輸入事件')throw Error('Raw 2.45 input boundary incorrect');
// Emulate native Home then repeated ArrowRight steps through the actual DOM
// reflection helper; each 0.1-percent input must survive a redraw.
let reflected=0;
for(let i=1;i<=1000;i++){reflected=api.progressValue((reflected+.1)/100);if(Math.abs(reflected-i/10)>1e-9)throw Error('Forward keyboard step lost at '+i);}
for(let i=999;i>=0;i--){reflected=api.progressValue((reflected-.1)/100);if(Math.abs(reflected-i/10)>1e-9)throw Error('Backward keyboard step lost at '+i);}
// A repeated browser asset must reuse the original API before allocating any
// observer or scheduling a DOM-ready handler. The pure Node renderer above
// simultaneously verifies that the document-free export route is unaffected.
const vm=require('vm');let remounts=0;const fakeDoc={};const existing={mount(doc){if(doc!==fakeDoc)throw Error('Singleton remounted wrong document');remounts++;}};
const browser={window:{EarthquakeExhibits:existing},document:fakeDoc};
vm.runInNewContext(fs.readFileSync('book/_static/exhibits.js','utf8'),browser);
vm.runInNewContext(fs.readFileSync('book/_static/exhibits.js','utf8'),browser);
if(browser.window.EarthquakeExhibits!==existing||remounts!==2)throw Error('Browser singleton replaced API');
const protocol=api.render('protocol',.6,{variant:'共同時鐘'}, {},390);
if(!protocol.body.includes('評估既有資料')||protocol.body.includes('事後選規則')||!protocol.caption.includes('樣本內或樣本外')||!protocol.caption.includes('子類')||!protocol.caption.includes('修訂 HORUS'))throw Error('Protocol retrospective semantics inconsistent');
const delayed=api.render('protocol',.6,{variant:'延遲驗證'}, {},390);
if(!delayed.caption.includes('模型及建模者')||!delayed.caption.includes('受控隔離'))throw Error('Delayed prospective isolation missing');
const partial=api.render('protocol',.6,{variant:'部分前瞻'}, {},390);
if(!partial.body.includes('已知')||!partial.body.includes('未知')||!partial.caption.includes('部分前瞻'))throw Error('Partial prospective distinction missing');
const hawkesBoundaries=[];
for(const [time,n] of [[2,0],[5,1],[5.8,2]]){
 const at=api.hawkesAt(time),r=api.render('hawkes',time/8,{}, {},760);
 const rasterCount=(r.body.match(/class="ex-event"/g)||[]).length;
 if(at.history.length!==n||rasterCount!==n||r.value!==at.total.toFixed(2)+' 顆／天')throw Error('Hawkes left-limit boundary mismatch at '+time);
 hawkesBoundaries.push({time,history_count:at.history.length,raster_count:rasterCount,total:at.total,passed:true});
}
const single=api.hawkesAt(3),singleSvg=api.render('hawkes',3/8,{}, {},760);
if(single.components.length!==1||Math.abs(single.background+single.components[0]-single.total)>1e-12)throw Error('Hawkes component sum mismatch');
const firstComponent=singleSvg.body.match(/<path d="M[0-9.]+ ([0-9.]+)[^"]*" class="ex-component ex-component-0"/);
const baseline=singleSvg.height-105,scale=(baseline-40)/4.3;
if(!firstComponent||Math.abs(Number(firstComponent[1])-(baseline-1.5*scale))>.011)throw Error('Hawkes component includes background in SVG');
const files=['book/_static/exhibits.js','book/_static/exhibits.css','book/_static/diagrams/sources/exhibits.json'];
const report={checked_at_utc:new Date().toISOString(),checks,total_states:checks.reduce((n,r)=>n+r.states,0),hawkes_semantics:{exact_boundaries:hawkesBoundaries,one_event_component_plus_background:single,component_svg_zero_baseline:true},protocol_semantics:{retrospective_in_or_out_of_sample:true,pseudo_is_subset:true,delayed_isolation:true,partial_known_unknown:true,revised_catalogue_limit:true},browser_singleton:{repeated_loads:2,existing_api_preserved:true,mount_calls:2,node_export_preserved:true},slider_precision:{step_percent:.1,forward_steps:1000,backward_steps:1000,passed:true},math_examples:['raw 2.45 enters input','Hawkes background .3 before event','fixed mixture 14','discrete upper tail .85','raw 4.95 enters target'],data:{cells:177,history_events:128,targets:25,windows:40},files:files.map(file=>({file,sha256:crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex')})),browser_status:'not run; main agent coordinates'};
fs.writeFileSync(path.resolve(__dirname,'../standalone/checks/exhibits_static.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({exhibits:checks.length,states:report.total_states,passed:true}));
