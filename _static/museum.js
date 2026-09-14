/* Navigation, reading progress and a bounded replay of actual target events. */
(() => {
  function start() {
    const shell=document.querySelector('.museum-shell');if(!shell)return;
    const slot=shell.querySelector('[data-museum-theme-slot]');
    const theme=document.querySelector('button.theme-switch-button');
    if(slot&&theme){slot.appendChild(theme);theme.setAttribute('aria-label','切換深淺主題');theme.title='切換深淺主題';}
    const sidebar=document.querySelector('.bd-sidebar-primary');
    const current=sidebar?.querySelector('a.current')||[...sidebar?.querySelectorAll('.bd-sidenav a')||[]].find(a=>new URL(a.href,location.href).pathname===location.pathname);
    if(current){
      current.setAttribute('aria-current','page');current.classList.add('current');
      const revealCurrent=()=>{if(matchMedia('(min-width:960px)').matches){const a=current.getBoundingClientRect(),b=sidebar.getBoundingClientRect();sidebar.scrollTop+=a.top-b.top-sidebar.clientHeight/3;}};
      requestAnimationFrame(revealCurrent);
      matchMedia('(min-width:960px)').addEventListener('change',revealCurrent);
    }
    const dialog=shell.querySelector('.museum-menu');
    matchMedia('(min-width:960px)').addEventListener('change',e=>{if(e.matches&&dialog.open)dialog.close();});
    shell.querySelector('[data-museum-menu]').addEventListener('click',()=>dialog.showModal());
    shell.querySelector('[data-museum-close]').addEventListener('click',()=>dialog.close());
    dialog.addEventListener('click',e=>{if(e.target===dialog){const r=dialog.getBoundingClientRect();if(e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom)dialog.close();}});
    shell.querySelector('.museum-search').addEventListener('click',event=>{const search=document.querySelector('#pst-search-dialog');if(search?.showModal){event.preventDefault();search.showModal();search.querySelector('input')?.focus();}});
    let pending=false;const progress=shell.querySelector('.museum-reading-progress');
    const update=()=>{pending=false;const maximum=document.documentElement.scrollHeight-innerHeight;progress.style.width=(maximum>0?Math.max(0,Math.min(100,scrollY/maximum*100)):0)+'%';};
    addEventListener('scroll',()=>{if(!pending){pending=true;requestAnimationFrame(update);}},{passive:true});addEventListener('resize',update);update();
    const source=document.querySelector('#museum-atlas-data'),svg=document.querySelector('[data-museum-atlas]');
    if(!source||!svg)return;
    const data=JSON.parse(source.textContent);if(!data.polygon?.length)return;
    const ns='http://www.w3.org/2000/svg',year=document.querySelector('[data-atlas-year]'),range=document.querySelector('[data-atlas-range]'),counter=document.querySelector('[data-atlas-count]'),button=document.querySelector('[data-atlas-play]');
    const points=data.polygon,xs=points.map(p=>p[0]),ys=points.map(p=>p[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    const scale=Math.min(450/(xmax-xmin),520/(ymax-ymin));const X=x=>35+(450-(xmax-xmin)*scale)/2+(x-xmin)*scale,Y=y=>35+(ymax-y)*scale;
    function make(tag,attrs){const el=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([key,value])=>el.setAttribute(key,value));return el;}
    svg.appendChild(make('title',{})).textContent='義大利實驗格網與目標地震回放';
    svg.appendChild(make('path',{d:points.map((p,i)=>(i?'L':'M')+X(p[0])+','+Y(p[1])).join(' ')+'Z',class:'atlas-collection'}));
    data.cells.forEach(c=>svg.appendChild(make('rect',{x:X(c[0]),y:Y(c[3]),width:(c[1]-c[0])*scale,height:(c[3]-c[2])*scale,class:'atlas-cell'})));
    const events=make('g',{});svg.appendChild(events);
    let timer=null,previous=2021;
    function draw(){const selected=+range.value;year.textContent=selected;events.replaceChildren();const visible=data.targets.filter(e=>e[2]<=selected);visible.forEach(e=>{const dot=make('circle',{cx:X(e[0]),cy:Y(e[1]),r:3+(e[3]-5)*3,class:'atlas-event'+(e[2]>previous?' atlas-new':'')});const title=make('title',{});title.textContent=`${e[4].slice(0,10)} · Mw ${e[3].toFixed(1)}`;dot.appendChild(title);events.appendChild(dot);});counter.textContent=visible.length+' 顆';previous=selected;}
    function pause(){if(timer)clearInterval(timer);timer=null;button.textContent='播放回放';button.setAttribute('aria-label','播放目標地震回放');}
    button.addEventListener('click',()=>{if(timer){pause();return;}if(+range.value>=2021)range.value=2012;previous=2011;draw();button.textContent='暫停';button.setAttribute('aria-label','暫停目標地震回放');timer=setInterval(()=>{if(+range.value>=2021){pause();return;}range.value=+range.value+1;draw();},800);});
    range.addEventListener('input',()=>{pause();draw();});new IntersectionObserver(entries=>{if(!entries[0].isIntersecting)pause();}).observe(svg);document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
    draw();
  }
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',start);else start();
})();
