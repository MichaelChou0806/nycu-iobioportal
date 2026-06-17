// =====================================================================
// analyses/immuneCorr.js —「Immune Correlation」
// ---------------------------------------------------------------------
// GOI 表現 × 免疫分數 相關（Spearman 預設 / Pearson）。只取 tumor。
// 篩選器（三件式）：方法快選(TIMER/CIBERSORT/.../EPIC) + 關鍵字搜尋(CD8/CAF…) + 已選清單。
// 通用維度選擇器（基因 / 癌種 / 細胞類型 三維）：
//   全單 → 散點圖 ; 一維多 → r-bar(+95% CI) ; 兩維多 → heatmap(可交換行列) ; 三維多 → 提示收一維
// r 的 p 與 CI 用 Fisher z。per-cancer 計算（不混癌種）。GOIs 與其他分析共享。
// =====================================================================

import { patientsInScope } from "../core/dimensions.js";
import { pearsonr, spearmanr, benjaminiHochberg, pStars } from "../core/stats.js";
import { corrScatterSVG, corrBarSVG, heatmapSVG } from "../core/plots.js";
import { getGOIs, setGOIs, onGOIsChanged, parseGenes } from "../core/gois.js";

const DIMS = ["gene", "cancer", "cell"];
// 顯示用：把 "<cell>_<METHOD>" 顯示成 "<cell> (METHOD)"
function formatCell(name) { const i = name.lastIndexOf("_"); return i < 0 ? name : name.slice(0, i) + " (" + name.slice(i + 1) + ")"; }
// 標題用：cell 拆成 [名稱, (方法)]（在「方法」前折行），方便標題分行
function cellTitleLines(rawCell) { const fc = formatCell(rawCell); const i = fc.lastIndexOf(" ("); return i < 0 ? [fc] : [fc.slice(0, i), fc.slice(i + 1)]; }

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .im-sec{margin-bottom:18px}
    .im-h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .im-sum{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;font-weight:600;cursor:pointer;padding:4px 0}
    details.im-sec>summary{margin-bottom:6px}
    .im-genes{width:100%;min-height:52px;padding:9px 11px;border:1px solid var(--line);border-radius:7px;font:13px/1.5 monospace;resize:vertical}
    .im-toolbar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;align-items:center}
    .im-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .im-mini.on{background:var(--accent);color:#fff}
    .im-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .im-chips{display:flex;flex-wrap:wrap;gap:6px}
    .im-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .im-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .im-chip.cell{background:#eef2ff;color:#3730a3;border-color:#c7d2fe;cursor:default}
    .im-chip.cell b{cursor:pointer;margin-left:4px;font-weight:700}
    .im-chip.dragging{opacity:.4}
    .im-search{width:260px;max-width:100%;padding:6px 10px;font-size:12px;border:1px solid var(--line);border-radius:7px}
    .im-results{display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:5px 14px;margin:8px 0;max-height:170px;overflow:auto;padding-right:6px}
    .im-srow{display:flex;align-items:flex-start;gap:6px;font-size:12px;color:#1f2733;cursor:pointer;line-height:1.4;min-width:0}
    .im-srow input{margin-top:2px;flex:0 0 auto}
    .im-selbox{display:flex;flex-wrap:wrap;gap:6px;margin-top:6px}
    .im-note{font-size:11.5px;color:var(--muted)}
    .im-ctl{font-size:11.5px;color:var(--muted)}
    .im-ctl select,.im-ctl input{padding:5px 8px;font-size:12px;margin-left:4px}
    .im-result{overflow-x:auto}
    .im-result svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .im-legend{font-size:11.5px;color:var(--muted);margin-top:6px}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

const LS = "tcga-tool:immune";
function loadState() { try { return JSON.parse(localStorage.getItem(LS)) || {}; } catch (e) { return {}; } }
function saveState(st) { try { localStorage.setItem(LS, JSON.stringify(st)); } catch (e) {} }
function moveInArray(arr, fromId, toId) {
  if (fromId == null || fromId === toId) return arr.slice();
  const a = arr.slice(); const fi = a.indexOf(fromId); if (fi < 0) return a;
  a.splice(fi, 1); let ti = a.indexOf(toId); if (ti < 0) ti = a.length; a.splice(ti, 0, fromId); return a;
}
// r 的 95% CI（Fisher z）
function corrCI(r, n) {
  if (n < 4 || Math.abs(r) >= 1) return { lo: NaN, hi: NaN };
  const z = 0.5 * Math.log((1 + r) / (1 - r)), se = 1 / Math.sqrt(n - 3);
  return { lo: Math.tanh(z - 1.96 * se), hi: Math.tanh(z + 1.96 * se) };
}
function fmtP(p) { return p < 0.001 ? "< 0.001" : "= " + p.toFixed(3); }

export const immuneCorr = {
  id: "immuneCorr",
  name: "Immune Correlation",

  async mount(container, { dataset }) {
    injectStyles();
    const availCancers = dataset.cancers.map(c => c.code);
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));
    const saved = loadState();
    let state = {
      cancers: (Array.isArray(saved.cancers) ? saved.cancers.filter(c => availCancers.includes(c)) : []).concat(availCancers.filter(c => !(saved.cancers || []).includes(c))),
      selectedCancers: Array.isArray(saved.selectedCancers) ? saved.selectedCancers.filter(c => availCancers.includes(c)) : [],
      selectedCellTypes: Array.isArray(saved.selectedCellTypes) ? saved.selectedCellTypes : [],
      corrMethod: saved.corrMethod || "spearman",
      scheme: saved.scheme || "rb",
      colorMax: saved.colorMax != null ? saved.colorMax : 0,
      swap: !!saved.swap,
      onedView: saved.onedView || "heatmap",
      barYMax: saved.barYMax != null ? saved.barYMax : 0,                // 0 = 自動縮放
      labelAngle: saved.labelAngle != null ? saved.labelAngle : 90,      // 預設垂直
    };
    let dragCancer = null, lastSVGName = "immune", immune = null, geneVals = null, byCancer = null;

    container.innerHTML = `
      <div class="card">
        <div class="im-sec">
          <h3 class="im-h3">Genes (shared with other analyses)</h3>
          <textarea id="im-genes" class="im-genes"></textarea>
        </div>

        <details class="im-sec" open>
          <summary class="im-sum">Cancers — drag to reorder · per-cancer correlation</summary>
          <div class="im-toolbar" style="margin-top:8px">
            <button class="im-mini" id="im-alpha">Sort A–Z</button>
            <button class="im-mini" id="im-byn">Sort by N</button>
            <button class="im-mini" id="im-def">Default</button>
            <span class="im-sep"></span>
            <button class="im-mini" id="im-call">Select all</button>
            <button class="im-mini" id="im-cnone">Select none</button>
          </div>
          <div id="im-cancers" class="im-chips"></div>
        </details>

        <details class="im-sec" open>
          <summary class="im-sum">Immune cell types — pick a whole method or search a keyword</summary>
          <div id="im-methods" class="im-toolbar" style="margin-top:8px"></div>
          <div class="im-toolbar">
            <input id="im-search" class="im-search" placeholder="search e.g. CD8, fibroblast, CAF, macrophage…">
            <button class="im-mini" id="im-clearcells">Clear selected</button>
            <span class="im-note" id="im-count">0 selected</span>
          </div>
          <div id="im-results" class="im-results"></div>
          <div id="im-selected" class="im-selbox"></div>
        </details>

        <div class="im-sec">
          <div class="im-toolbar">
            <label class="im-ctl">Method
              <select id="im-method"><option value="spearman">Spearman</option><option value="pearson">Pearson</option></select></label>
            <label class="im-ctl" id="im-oned-wrap">1-D view
              <select id="im-oned"><option value="heatmap">Heatmap</option><option value="bar">Bar</option></select></label>
            <label class="im-ctl" id="im-scheme-wrap">Colors
              <select id="im-scheme"><option value="rb">Red–Blue</option><option value="rg">Red–Green</option></select></label>
            <label class="im-ctl" id="im-cmax-wrap">Scale max<input type="number" id="im-cmax" min="0" max="1" step="0.1" style="width:60px" placeholder="1"></label>
            <label class="im-ctl" id="im-ymax-wrap">Y max<input type="number" id="im-ymax" min="0.05" max="1" step="0.05" style="width:60px" placeholder="auto"></label>
            <label class="im-ctl" id="im-angle-wrap">Label °<input type="number" id="im-angle" min="0" max="90" step="5" style="width:54px"></label>
            <button class="im-mini" id="im-swap" title="swap heatmap rows/cols">Swap rows/cols</button>
            <button id="im-run">Run</button>
            <button class="im-mini" id="im-svg" style="display:none">Download SVG</button>
            <button class="im-mini" id="im-png" style="display:none">Download PNG</button>
          </div>
          <div id="im-status" class="status">Loading immune data…</div>
          <div id="im-result" class="im-result"></div>
        </div>
        <div class="im-note">Tumor only · per-cancer (no cross-cancer pooling). r-bar error bars = 95% CI (Fisher z). Heatmap/bar p are FDR-corrected.</div>
      </div>`;

    const $ = s => container.querySelector(s);
    const genesEl = $("#im-genes"), chipBox = $("#im-cancers"), statusEl = $("#im-status"), resultEl = $("#im-result");
    const methodBox = $("#im-methods"), searchInput = $("#im-search"), searchResults = $("#im-results"), selectedBox = $("#im-selected"), countEl = $("#im-count");
    const methodSel = $("#im-method"), schemeSel = $("#im-scheme"), cmaxEl = $("#im-cmax"), ymaxEl = $("#im-ymax"), angleEl = $("#im-angle");

    function commit() { saveState(state); }
    genesEl.value = getGOIs().join(", ");
    genesEl.addEventListener("input", () => setGOIs(parseGenes(genesEl.value)));
    onGOIsChanged(list => { if (document.activeElement !== genesEl) genesEl.value = list.join(", "); updateControls(); });
    methodSel.value = state.corrMethod; schemeSel.value = state.scheme; $("#im-oned").value = state.onedView; if (state.colorMax > 0) cmaxEl.value = state.colorMax;
    if (state.barYMax > 0) ymaxEl.value = state.barYMax; angleEl.value = state.labelAngle;

    // 控件可見性（依當前選到幾個維度為「多」）
    function show(sel, on) { const el = $(sel); if (el) el.style.display = on ? "" : "none"; }
    function updateControls() {
      const nG = getGOIs().length, nC = state.selectedCancers.length, nT = state.selectedCellTypes.length;
      const nMulti = (nG > 1 ? 1 : 0) + (nC > 1 ? 1 : 0) + (nT > 1 ? 1 : 0);
      show("#im-oned-wrap", nMulti === 1);                                       // 一維多才需選 bar/heatmap
      show("#im-scheme-wrap", nMulti >= 1);                                      // 散點不需配色
      show("#im-cmax-wrap", nMulti === 2 || (nMulti === 1 && state.onedView === "heatmap"));
      show("#im-ymax-wrap", nMulti === 1 && state.onedView === "bar");           // Y 軸縮放只在 bar
      show("#im-angle-wrap", (nMulti === 1 && state.onedView === "bar") || nMulti === 2);  // 角度：bar 與 2D heatmap（有旋轉長標籤）
      show("#im-swap", nMulti === 2);                                            // 只有兩維 heatmap 能交換
    }

    // ---- 癌種 chips ----
    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "im-chip" + (state.selectedCancers.includes(code) ? " on" : "");
        chip.textContent = code; chip.draggable = true;
        chip.addEventListener("click", () => {
          if (state.selectedCancers.includes(code)) state.selectedCancers = state.selectedCancers.filter(c => c !== code);
          else state.selectedCancers.push(code);
          commit(); renderCancers();
        });
        chip.addEventListener("dragstart", () => { dragCancer = code; chip.classList.add("dragging"); });
        chip.addEventListener("dragend", () => chip.classList.remove("dragging"));
        chip.addEventListener("dragover", e => e.preventDefault());
        chip.addEventListener("drop", e => { e.preventDefault(); state.cancers = moveInArray(state.cancers, dragCancer, code); commit(); renderCancers(); });
        chipBox.appendChild(chip);
      });
      updateControls();
    }
    $("#im-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#im-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#im-def").addEventListener("click", () => { state.cancers = availCancers.slice(); commit(); renderCancers(); });
    $("#im-call").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); });
    $("#im-cnone").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); });

    // ---- 細胞類型篩選器（需 immune 載入後才有清單）----
    function renderMethods() {
      methodBox.innerHTML = "";
      if (!immune) return;
      Object.entries(immune.methods).forEach(([m, cts]) => {
        const allOn = cts.every(c => state.selectedCellTypes.includes(c));
        const btn = document.createElement("button");
        btn.className = "im-mini" + (allOn ? " on" : "");
        btn.textContent = `${m} (${cts.length})`;
        btn.addEventListener("click", () => {
          if (allOn) state.selectedCellTypes = state.selectedCellTypes.filter(c => !cts.includes(c));
          else cts.forEach(c => { if (!state.selectedCellTypes.includes(c)) state.selectedCellTypes.push(c); });
          commit(); renderMethods(); renderSelected(); renderSearch();
        });
        methodBox.appendChild(btn);
      });
    }
    function renderSearch() {
      searchResults.innerHTML = "";
      if (!immune) return;
      const kw = searchInput.value.trim().toLowerCase();
      if (!kw) { searchResults.innerHTML = `<span class="im-note">Type a keyword (CD8, fibroblast, CAF, macrophage, NK…) to find cell types across all methods.</span>`; return; }
      const matches = immune.cellTypes.filter(c => c.toLowerCase().includes(kw)).slice(0, 80);
      if (!matches.length) { searchResults.innerHTML = `<span class="im-note">No cell type contains “${kw}”.</span>`; return; }
      matches.forEach(c => {
        const on = state.selectedCellTypes.includes(c);
        const row = document.createElement("label"); row.className = "im-srow";
        const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = on;
        const sp = document.createElement("span"); sp.textContent = formatCell(c);
        cb.addEventListener("change", () => {
          if (cb.checked) { if (!state.selectedCellTypes.includes(c)) state.selectedCellTypes.push(c); }
          else state.selectedCellTypes = state.selectedCellTypes.filter(x => x !== c);
          commit(); renderSelected(); renderMethods();
        });
        row.appendChild(cb); row.appendChild(sp); searchResults.appendChild(row);
      });
    }
    function renderSelected() {
      selectedBox.innerHTML = "";
      countEl.textContent = `${state.selectedCellTypes.length} selected`;
      state.selectedCellTypes.forEach(c => {
        const chip = document.createElement("span"); chip.className = "im-chip cell";
        const t = document.createTextNode(formatCell(c) + " ");
        const x = document.createElement("b"); x.textContent = "×";
        x.addEventListener("click", () => { state.selectedCellTypes = state.selectedCellTypes.filter(v => v !== c); commit(); renderSelected(); renderMethods(); renderSearch(); });
        chip.appendChild(t); chip.appendChild(x); selectedBox.appendChild(chip);
      });
      updateControls();
    }
    searchInput.addEventListener("input", renderSearch);
    $("#im-clearcells").addEventListener("click", () => { state.selectedCellTypes = []; commit(); renderSelected(); renderMethods(); renderSearch(); });

    methodSel.addEventListener("change", () => { state.corrMethod = methodSel.value; commit(); });
    schemeSel.addEventListener("change", () => { state.scheme = schemeSel.value; commit(); updateControls(); });
    $("#im-oned").addEventListener("change", () => { state.onedView = $("#im-oned").value; commit(); updateControls(); });
    cmaxEl.addEventListener("change", () => { const v = Number(cmaxEl.value); state.colorMax = (cmaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : Math.min(1, v); commit(); });
    ymaxEl.addEventListener("change", () => { const v = Number(ymaxEl.value); state.barYMax = (ymaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : Math.min(1, v); commit(); if (resultEl.querySelector("svg")) run(); });
    angleEl.addEventListener("change", () => { let v = Number(angleEl.value); if (!isFinite(v)) v = 90; v = Math.max(0, Math.min(90, v)); state.labelAngle = v; angleEl.value = v; commit(); if (resultEl.querySelector("svg")) run(); });
    $("#im-swap").addEventListener("click", () => { state.swap = !state.swap; commit(); if (resultEl.querySelector("svg")) run(); });

    // ---- 載入免疫資料 ----
    try {
      immune = await dataset.loadImmune();
      statusEl.className = "status"; statusEl.textContent = `Immune data ready: ${immune.cellTypes.length} cell types.`;
      renderMethods(); renderSearch(); renderSelected();
    } catch (e) {
      statusEl.className = "status err"; statusEl.textContent = "Failed to load immune data (CORS / not uploaded?): " + e.message;
    }

    // ---- 計算單一組合的相關 ----
    function corrOf(geneObj, cancer, cellType) {
      const exprArr = geneVals.get(geneObj.rec.gene_id), immuneArr = immune.get(cellType);
      const ps = byCancer[cancer] || [];
      const xs = [], ys = [];
      for (const p of ps) { const e = exprArr[p.idx], m = immuneArr ? immuneArr[p.idx] : NaN; if (isFinite(e) && isFinite(m)) { xs.push(m); ys.push(e); } }
      const n = xs.length;
      if (n < 4) return { n, state: "nodata", tip: `${geneObj.label} / ${cancer} / ${cellType}: n=${n}` };
      const res = state.corrMethod === "pearson" ? pearsonr(xs, ys) : spearmanr(xs, ys);
      const ci = corrCI(res.r, n);
      const weak = n < 20;
      return { r: res.r, p: res.p, ciLow: ci.lo, ciHigh: ci.hi, n, xs, ys, state: weak ? "weak" : "ok", tip: `${geneObj.label} / ${cancer} / ${cellType}: r=${res.r.toFixed(3)}, p=${res.p.toPrecision(2)}, n=${n}` };
    }
    const methodName = () => state.corrMethod === "pearson" ? "Pearson" : "Spearman";

    // ---- RUN ----
    $("#im-run").addEventListener("click", run);
    $("#im-svg").addEventListener("click", downloadSVG);
    $("#im-png").addEventListener("click", downloadPNG);

    async function run() {
      resultEl.innerHTML = ""; $("#im-svg").style.display = "none"; $("#im-png").style.display = "none";
      if (!immune) { statusEl.className = "status err"; statusEl.textContent = "Immune data not loaded."; return; }
      const cancers = state.cancers.filter(c => state.selectedCancers.includes(c));
      const cellTypes = state.selectedCellTypes.filter(c => immune.has(c));
      if (!cancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      if (!cellTypes.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one immune cell type."; return; }
      const recs = []; const unknown = [];
      getGOIs().forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) unknown.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (!recs.length) { statusEl.className = "status err"; statusEl.textContent = "No recognized genes." + (unknown.length ? ` Unrecognized: ${unknown.join(", ")}` : ""); return; }

      statusEl.className = "status"; statusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      try { geneVals = new Map(); await Promise.all(recs.map(async x => geneVals.set(x.rec.gene_id, await dataset.getGeneValues(x.rec)))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files: " + e.message; return; }

      statusEl.textContent = "Computing…";
      const patients = patientsInScope(dataset, cancers); // tumor only
      byCancer = {}; cancers.forEach(c => byCancer[c] = patients.filter(p => p.cancer === c));

      const itemsOf = d => d === "gene" ? recs.map(r => ({ val: r, label: r.label })) : d === "cancer" ? cancers.map(c => ({ val: c, label: c })) : cellTypes.map(c => ({ val: c, label: formatCell(c) }));
      const single = { gene: recs[0], cancer: cancers[0], cell: cellTypes[0] };
      const multiKeys = DIMS.filter(d => itemsOf(d).length > 1);

      if (multiKeys.length === 0) drawScatter(recs[0], cancers[0], cellTypes[0]);
      else if (multiKeys.length === 1) { if (state.onedView === "bar") drawBar(multiKeys[0], itemsOf, single); else drawOneDimHeatmap(multiKeys[0], itemsOf, single); }
      else if (multiKeys.length === 2) drawHeatmap(multiKeys, itemsOf, single);
      else { statusEl.className = "status err"; statusEl.textContent = "All three of genes / cancers / cell types have multiple selections. Please reduce one of them to a single item (e.g. one gene, or one cancer) so it can be plotted."; return; }

      const note = document.createElement("div"); note.className = "im-legend";
      note.textContent = `${methodName()} · tumor only · per-cancer${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      $("#im-svg").style.display = ""; $("#im-png").style.display = "";
      statusEl.textContent = "Done.";
    }

    function drawScatter(geneObj, cancer, cellType) {
      lastSVGName = `corr_${geneObj.label}_${cancer}`;
      const c = corrOf(geneObj, cancer, cellType);
      if (c.state === "nodata") { resultEl.innerHTML = `<div class="status err">Too few paired samples (n=${c.n}).</div>`; return; }
      const points = c.xs.map((x, i) => ({ x, y: c.ys[i] }));
      resultEl.innerHTML = corrScatterSVG(points, { xLabel: formatCell(cellType), yLabel: `${geneObj.label} expression`, r: c.r, pText: fmtP(c.p), n: c.n, title: `${geneObj.label} vs ${formatCell(cellType)}`, subtitle: cancer });
    }

    function drawBar(md, itemsOf, single) {
      lastSVGName = `corr_bar_${md}`;
      const items = itemsOf(md);
      const computed = items.map(it => { const dv = { ...single }; dv[md] = it.val; return { label: it.label, c: corrOf(dv.gene, dv.cancer, dv.cell) }; });
      const pv = computed.filter(x => x.c.state !== "nodata").map(x => x.c.p);
      const q = benjaminiHochberg(pv); let qi = 0;
      computed.forEach(x => { if (x.c.state !== "nodata") { const qq = q[qi++]; const st = pStars(qq); x.c.stars = st === "ns" ? "" : st; x.c.tip += `, q=${qq.toPrecision(2)}`; } });
      const barItems = computed.map(x => x.c.state === "nodata" ? { label: x.label, state: "nodata" } : { label: x.label, r: x.c.r, ciLow: x.c.ciLow, ciHigh: x.c.ciHigh, stars: x.c.stars, state: x.c.state, tip: x.c.tip });
      const capLines = [];                                            // Y 軸已是「Spearman r」，標題不重複 method
      ["cell", "gene", "cancer"].filter(d => d !== md).forEach(d => {
        if (d === "cell") capLines.push(...cellTitleLines(single.cell));
        else if (d === "gene") capLines.push(single.gene.label);
        else capLines.push(single.cancer);
      });
      resultEl.innerHTML = corrBarSVG(barItems, { caption: capLines.join("\n"), yLabel: `${methodName()} r`, scheme: state.scheme, yMax: state.barYMax, labelAngle: state.labelAngle });
    }

    function drawOneDimHeatmap(md, itemsOf, single) {
      lastSVGName = `corr_1d_${md}`;
      const rowItems = itemsOf(md);
      // X 列標：優先固定的癌種，其次基因，其次細胞（避免無意義的 "r"）
      const xDim = (md !== "cancer") ? "cancer" : (md !== "gene" ? "gene" : "cell");
      const xLabel = xDim === "cancer" ? single.cancer : xDim === "gene" ? single.gene.label : formatCell(single.cell);
      // 標題：method correlation + 固定的描述維度（非 md、非 X）；不重複 Y 軸已有的 md 維度
      const titleDim = DIMS.find(d => d !== md && d !== xDim);
      let capLines = [`${methodName()} correlation`];
      if (titleDim === "cell") capLines = capLines.concat(cellTitleLines(single.cell));
      else if (titleDim === "gene") capLines.push(single.gene.label);
      else capLines.push(single.cancer);
      const cells = []; const pflat = [];
      rowItems.forEach((it, r) => {
        const dv = { ...single }; dv[md] = it.val;
        const c = corrOf(dv.gene, dv.cancer, dv.cell);
        if (c.state === "nodata") cells.push([{ state: "nodata", tip: c.tip }]);
        else { cells.push([{ value: c.r, state: c.state, stars: "", tip: c.tip }]); pflat.push({ r, p: c.p }); }
      });
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = cells[x.r][0]; const qq = q[k]; const st = pStars(qq); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${qq.toPrecision(2)}`; });
      const colorMax = state.colorMax > 0 ? state.colorMax : 1;
      resultEl.innerHTML = heatmapSVG(rowItems.map(x => x.label), [xLabel], cells, { scheme: state.scheme, colorMax, legendLabel: "r", caption: capLines.join("\n") });
    }

    function drawHeatmap(multiKeys, itemsOf, single) {
      lastSVGName = "corr_heatmap";
      let [rowDim, colDim] = multiKeys; if (state.swap) [rowDim, colDim] = [colDim, rowDim];
      const fixedDim = DIMS.find(d => !multiKeys.includes(d));
      const rowItems = itemsOf(rowDim), colItems = itemsOf(colDim);
      const cells = []; const pflat = [];
      rowItems.forEach((ri, r) => {
        const rowCells = [];
        colItems.forEach((ci, cc) => {
          const dv = { ...single }; dv[rowDim] = ri.val; dv[colDim] = ci.val;
          const c = corrOf(dv.gene, dv.cancer, dv.cell);
          if (c.state === "nodata") rowCells.push({ state: "nodata", tip: c.tip });
          else { rowCells.push({ value: c.r, state: c.state, stars: "", tip: c.tip }); pflat.push({ r, cc, p: c.p }); }
        });
        cells.push(rowCells);
      });
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = cells[x.r][x.cc]; const qq = q[k]; const st = pStars(qq); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${qq.toPrecision(2)}`; });
      const colorMax = state.colorMax > 0 ? state.colorMax : 1;
      const fixedLabel = fixedDim === "gene" ? single.gene.label : fixedDim === "cancer" ? single.cancer : formatCell(single.cell);
      resultEl.innerHTML = heatmapSVG(rowItems.map(x => x.label), colItems.map(x => x.label), cells, { scheme: state.scheme, colorMax, legendLabel: "r", caption: `${methodName()} correlation · ${fixedLabel}`, labelAngle: state.labelAngle });
    }

    // ---- 下載 ----
    function firstSVG() { return resultEl.querySelector("svg"); }
    function dl(blob, name) { const u = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = u; a.download = name; a.click(); URL.revokeObjectURL(u); }
    function downloadSVG() { const svg = firstSVG(); if (!svg) return; dl(new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml;charset=utf-8" }), lastSVGName + ".svg"); }
    function downloadPNG() {
      const svg = firstSVG(); if (!svg) return;
      const vb = svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal : null;
      const W = vb && vb.width ? vb.width : svg.clientWidth, H = vb && vb.height ? vb.height : svg.clientHeight;
      const scale = Math.max(2, Math.ceil(2400 / W));
      const xml = new XMLSerializer().serializeToString(svg);
      const img = new Image();
      img.onload = () => { const c = document.createElement("canvas"); c.width = Math.round(W * scale); c.height = Math.round(H * scale); const ctx = c.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, c.width, c.height); ctx.drawImage(img, 0, 0, c.width, c.height); c.toBlob(b => dl(b, lastSVGName + ".png"), "image/png"); };
      img.onerror = () => alert("PNG export failed (try SVG).");
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(xml);
    }

    renderCancers();
  },
};
