/* Keep saved Plotly figures in the same light/dark palette as the reading page. */
(() => {
  const applied = new WeakMap();
  const waitingMaps = new WeakSet();
  let scheduled = false;
  function updateFigures() {
    scheduled = false;
    if (!window.Plotly) return;
    const dark = document.documentElement.dataset.theme === 'dark';
    const key = dark ? 'dark' : 'light';
    const style=getComputedStyle(document.documentElement);
    const ink=style.getPropertyValue('--tc-ink').trim();
    const paper=style.getPropertyValue('--tc-paper').trim();
    const rule=style.getPropertyValue('--tc-rule').trim();
    document.querySelectorAll('.js-plotly-plot').forEach(fig => {
      if (!fig._fullLayout || applied.get(fig) === key) return;
      // MapLibre initializes its style asynchronously after Plotly's layout exists.
      // Relayout before its first idle event can access a style still being parsed.
      const mapKeys = Object.keys(fig._fullLayout).filter(name => /^map(box)?\d*$/.test(name));
      const maps = mapKeys.map(name => fig._fullLayout[name]._subplot?.map).filter(Boolean);
      if (maps.length < mapKeys.length) {
        setTimeout(schedule, 50);
        return;
      }
      const pending = maps.filter(map => !map.isStyleLoaded());
      if (pending.length) {
        pending.forEach(map => {
          if (waitingMaps.has(map)) return;
          waitingMaps.add(map);
          map.once('idle', () => { waitingMaps.delete(map); schedule(); });
        });
        return;
      }
      applied.set(fig, key);
      const layout = {
        paper_bgcolor: paper, plot_bgcolor: paper,
        'font.color': ink, 'legend.font.color': ink,
        'legend.bgcolor': 'rgba(0,0,0,0)',
      };
      Object.keys(fig._fullLayout).forEach(name => {
        if (/^[xy]axis\d*$/.test(name)) {
          layout[`${name}.color`] = ink;
          layout[`${name}.gridcolor`] = rule;
          layout[`${name}.zerolinecolor`] = rule;
        }
        if (/^geo\d*$/.test(name)) {
          layout[`${name}.bgcolor`] = paper;
          layout[`${name}.landcolor`] = dark ? '#2a313c' : '#eef2f7';
          layout[`${name}.coastlinecolor`] = dark ? '#b6c0cd' : '#56657a';
          layout[`${name}.countrycolor`] = dark ? '#b6c0cd' : '#56657a';
        }
      });
      Promise.resolve(Plotly.relayout(fig, layout)).catch(error => {
        console.error('Teaching figure theme failed', error);
      });
    });
  }
  function schedule() {
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(updateFigures);
    }
  }
  new MutationObserver(schedule).observe(document.documentElement,
    {attributes: true, attributeFilter: ['data-theme']});
  const start = () => {
    new MutationObserver(schedule).observe(document.body, {childList: true, subtree: true});
    schedule();
  };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start);
  else start();
  window.addEventListener('load', schedule);
})();
