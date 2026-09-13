/* Build static-first exhibit fragments using the same pure renderer as browser. */
const fs=require('fs'),path=require('path'),crypto=require('crypto');
const root=path.resolve(__dirname,'..'),api=require('../../exhibits.js');
const metadata=JSON.parse(fs.readFileSync(path.join(__dirname,'exhibits.json'),'utf8'));
const e=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const labels={S:'收集區 S',R:'177 個測試格',events:'歷史震央',contours:'等值圈',sources:'來源震央'};
const sourceLabels={real:'真實資料',synthetic:'固定合成示意',process:'實驗流程'};
const css=fs.readFileSync(path.resolve(root,'../exhibits.css'),'utf8');
const jscript=fs.readFileSync(path.resolve(root,'../exhibits.js'),'utf8');
const receipts=[];
for(const m of metadata){
 const def=api.kinds[m.kind],state={title:m.title,variant:def.select?def.select[0]:undefined},r=api.svgHTML(m.kind,1,state,m.data||{},760,m.slug+'-static');
 const layers=def.layers?`<div class="ex-layers">${def.layers.map(k=>`<label><input type="checkbox" data-layer="${k}" checked>${labels[k]}</label>`).join('')}</div>`:'';
 const select=def.select?`<label><span class="ex-sr">比較方式</span><select data-control="compare" aria-label="比較方式">${def.select.map(v=>`<option>${e(v)}</option>`).join('')}</select></label>`:'';
 const controls=`<div class="ex-controls">${def.play?'<button type="button" data-action="play" aria-pressed="false">播放導覽</button>':''}<button type="button" data-action="reset">從頭看</button>${select}${layers}${m.kind==='magnitude'?'<button type="button" data-action="boundary" data-value="2.45">看 2.45</button><button type="button" data-action="boundary" data-value="4.95">看 4.95</button>':''}<label class="ex-scrub">${m.control_label}<input type="range" data-control="progress" min="0" max="100" step="0.1" value="${api.progressValue(def.start)}" aria-label="${e(m.control_label)}"></label></div>`;
 const frag=`<section class="quake-exhibit" data-exhibit="${m.slug}" data-kind="${m.kind}" data-config="${e(JSON.stringify(m))}" aria-label="${e(m.title)}"><header class="ex-header"><div><p class="ex-kicker">${sourceLabels[m.source]} · ${m.slug.slice(1,3)}</p><h3 class="ex-title">${e(m.title)}</h3><p class="ex-deck">${e(m.deck)}</p></div><output class="ex-readout" aria-live="off">${e(r.value)}</output></header><div class="ex-stage">${r.svg}</div>${controls}<p class="ex-caption">${e(r.caption)}</p><span class="ex-status" role="status" aria-live="polite"></span></section>\n`;
 fs.writeFileSync(path.join(root,m.slug+'.html'),frag);
 const html=`<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>${e(m.title)}</title><style>:root{--tc-paper:#fcfcf9;--tc-paper-2:#eef3f2;--tc-ink:#203341;--tc-muted:#536773;--tc-rule:#ccd9d7;--tc-accent:#267e88;--tc-quake:#dc685c}html[data-theme=dark]{--tc-paper:#14252e;--tc-paper-2:#1c323c;--tc-ink:#e4edec;--tc-muted:#b5cbc9;--tc-rule:#3b555d;--tc-accent:#6ac4c3;--tc-quake:#f48b79}body{margin:0;padding:30px 20px;background:var(--tc-paper);color:var(--tc-ink)}main{max-width:1050px;margin:auto}${css}</style></head><body><main>${frag}</main><script>${jscript}</script></body></html>`;
 fs.writeFileSync(path.join(root,'standalone',m.slug+'.html'),html);
 receipts.push({slug:m.slug,kind:m.kind,fragment_sha256:crypto.createHash('sha256').update(frag).digest('hex'),standalone_sha256:crypto.createHash('sha256').update(html).digest('hex')});
}
fs.writeFileSync(path.join(root,'standalone/checks/exhibits_manifest.json'),JSON.stringify({count:metadata.length,generated_at_utc:new Date().toISOString(),exhibits:receipts,browser_status:'pending; main agent coordinates'},null,2)+'\n');
console.log(`Generated ${metadata.length} exhibits and offline standalone previews`);
