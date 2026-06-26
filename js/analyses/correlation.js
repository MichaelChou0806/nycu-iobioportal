// =====================================================================
// analyses/correlation.js —「Gene Correlation」
// ---------------------------------------------------------------------
// GOI 之間兩兩表現相關（Spearman 預設 / Pearson）。只取 tumor、per-cancer（不混癌種）。
// 路由（GOI 數 × 癌種數）：
//   2 基因 × 1 癌種 → scatter（回歸 + 95% 信賴帶 + 上/右邊緣分布）
//   2 基因 × 多癌種 → 1-D：r across cancers（heatmap / bar）
//   >2 基因 × 1 癌種 → 上三角相關矩陣（點格子 → 該對 scatter）
//   >2 基因 × 多癌種 → 提示收一軸（高維無法一張圖呈現，同 immuneCorr）
// miRNA（偽基因）自動可入：Spearman rank-based 讓 RPM×TPM 混合也合理。GOIs 與其他分析共享。
// =====================================================================

import { patientsInScope } from "../core/dimensions.js";
import { pearsonr, spearmanr, benjaminiHochberg, pStars } from "../core/stats.js";
import { corrScatterSVG, corrBarSVG, heatmapSVG, corrMatrixSVG } from "../core/plots.js";
import { getGOIs, setGOIs, onGOIsChanged, parseGenes } from "../core/gois.js";

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .cr-sec{margin-bottom:18px}
    .cr-h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .cr-sum{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;font-weight:600;cursor:pointer;padding:4px 0}
    details.cr-sec>summary{margin-bottom:6px}
    .cr-genes{width:100%;min-height:52px;padding:9px 11px;border:1px solid var(--line);border-radius:7px;font:13px/1.5 monospace;resize:vertical}
    .cr-toolbar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;align-items:center}
    .cr-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .cr-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .cr-chips{display:flex;flex-wrap:wrap;gap:6px}
    .cr-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .cr-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .cr-chip.dragging{opacity:.4}
    .cr-ctl{font-size:11.5px;color:var(--muted)}
    .cr-ctl select,.cr-ctl input{padding:5px 8px;font-size:12px;margin-left:4px}
    .cr-note{font-size:11.5px;color:var(--muted)}
    .cr-result{overflow-x:auto}
    .cr-result svg,.cr-pair svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .cr-pair{overflow-x:auto}
    .cr-legend{font-size:11.5px;color:var(--muted);margin-top:6px;line-height:1.5}
    .cr-views{display:flex;flex-wrap:wrap;gap:16px;align-items:flex-start}
    .cr-views > div{flex:1 1 360px;min-width:0}
    .cr-views > div:empty{display:none}
    .cr-pctl{display:flex;flex-wrap:wrap;gap:6px 10px;align-items:center;font-size:11.5px;color:var(--muted);margin-top:10px}
    .cr-pctl .cr-ptitle{flex-basis:100%;font-weight:600;color:#1f2733;font-size:12.5px}
    .cr-pn{width:62px;min-width:0;padding:4px 6px;font-size:12px;border:1px solid var(--line);border-radius:6px}
    .cr-popacity{width:90px;min-width:0;vertical-align:middle}
    .cr-plabel{padding:3px 6px;font-size:12px;border:1px solid var(--line);border-radius:6px;min-width:0}
    .cr-pswatches{display:flex;flex-wrap:wrap;gap:4px;margin:8px 0}
    .cr-sw{width:18px;height:18px;border-radius:4px;cursor:pointer;border:1px solid rgba(0,0,0,.12)}
    .cr-sw.on{outline:2px solid #1f2733;outline-offset:1px}
    .cr-pplot{overflow-x:auto}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

const LS = "tcga-tool:corr";
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
const unitOf = rec => rec.assay === "mirna_rpm" ? "RPM" : "TPM";
const csvq = s => /[",\n]/.test(String(s)) ? '"' + String(s).replace(/"/g, '""') + '"' : String(s);
// 散點填充色的基本色盤（24 色，外框由 plots.js 自動取同色加深）
const PALETTE = ["#3b82f6", "#2563eb", "#06b6d4", "#0891b2", "#14b8a6", "#10b981", "#22c55e", "#16a34a", "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444", "#dc2626", "#e11d48", "#ec4899", "#db2777", "#d946ef", "#a855f7", "#8b5cf6", "#6366f1", "#64748b", "#475569", "#1f2733"];

export const correlation = {
  id: "correlation",
  name: "Gene Correlation",

  async mount(container, { dataset }) {
    injectStyles();
    const availCancers = dataset.cancers.map(c => c.code);
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));
    const saved = loadState();
    let state = {
      cancers: (Array.isArray(saved.cancers) ? saved.cancers.filter(c => availCancers.includes(c)) : []).concat(availCancers.filter(c => !(saved.cancers || []).includes(c))),
      selectedCancers: Array.isArray(saved.selectedCancers) ? saved.selectedCancers.filter(c => availCancers.includes(c)) : (availCancers.includes("HNSC") ? ["HNSC"] : availCancers.slice(0, 1)),
      corrMethod: saved.corrMethod === "pearson" ? "pearson" : "spearman",
      scheme: saved.scheme === "rg" ? "rg" : "rb",
      colorMax: saved.colorMax != null ? saved.colorMax : 0,   // 0 = auto(=1)
      onedView: saved.onedView === "bar" ? "bar" : "heatmap",
      swap: !!saved.swap,
    };
    let dragCancer = null, lastSVGName = "corr", geneVals = null, byCancer = null, lastRecs = null, lastResult = null, lastScatterColor = "#3b82f6", lastScatterOpacity = 0.45, lastLabelPos = "tr";

    container.innerHTML = `
      <div class="card">
        <div class="cr-sec">
          <h3 class="cr-h3">Genes (shared with other analyses)</h3>
          <textarea id="cr-genes" class="cr-genes" placeholder="e.g. KDELR1, KDELR2, hsa-miR-21-5p …"></textarea>
        </div>

        <details class="cr-sec" open>
          <summary class="cr-sum">Cancers — drag to reorder · per-cancer correlation</summary>
          <div class="cr-toolbar" style="margin-top:8px">
            <button class="cr-mini" id="cr-alpha">Sort A–Z</button>
            <button class="cr-mini" id="cr-byn">Sort by N</button>
            <button class="cr-mini" id="cr-def">Default</button>
            <span class="cr-sep"></span>
            <button class="cr-mini" id="cr-call">Select all</button>
            <button class="cr-mini" id="cr-cnone">Select none</button>
          </div>
          <div id="cr-cancers" class="cr-chips"></div>
        </details>

        <div class="cr-sec">
          <div class="cr-toolbar">
            <label class="cr-ctl">Method<select id="cr-method"><option value="spearman">Spearman</option><option value="pearson">Pearson</option></select></label>
            <label class="cr-ctl" id="cr-oned-wrap">Multi-cancer<select id="cr-oned"><option value="heatmap">Heatmap</option><option value="bar">Bar</option></select></label>
            <label class="cr-ctl" id="cr-scheme-wrap">Colors<select id="cr-scheme"><option value="rb">Red–Blue</option><option value="rg">Red–Green</option></select></label>
            <label class="cr-ctl" id="cr-cmax-wrap">Color max<input type="number" id="cr-cmax" min="0.05" max="1" step="0.05" style="width:60px" placeholder="1"></label>
            <button class="cr-mini" id="cr-swap" title="swap rows/cols" style="display:none">Swap</button>
            <button id="cr-run">Run</button>
            <button class="cr-mini" id="cr-csv" style="display:none">CSV</button>
            <button class="cr-mini" id="cr-svg" style="display:none">SVG</button>
            <button class="cr-mini" id="cr-png" style="display:none">PNG</button>
          </div>
          <div id="cr-status" class="status"></div>
          <div class="cr-views">
            <div id="cr-result" class="cr-result"></div>
            <div id="cr-pair" class="cr-pair"></div>
          </div>
        </div>
        <div class="cr-note">Tumor only · per-cancer (no cross-cancer pooling). Pairwise-complete (drops samples missing either gene). Matrix shows the upper triangle — click a cell to see that pair's scatter. Heatmap/bar p are FDR-corrected.</div>
      </div>`;

    const $ = s => container.querySelector(s);
    const genesEl = $("#cr-genes"), chipBox = $("#cr-cancers"), statusEl = $("#cr-status"), resultEl = $("#cr-result"), pairEl = $("#cr-pair");
    const methodSel = $("#cr-method"), schemeSel = $("#cr-scheme"), cmaxEl = $("#cr-cmax"), onedSel = $("#cr-oned");

    function commit() { saveState(state); }
    genesEl.value = getGOIs().join(", ");
    genesEl.addEventListener("input", () => setGOIs(parseGenes(genesEl.value)));
    onGOIsChanged(list => { if (document.activeElement !== genesEl) genesEl.value = list.join(", "); updateControls(); });
    methodSel.value = state.corrMethod; schemeSel.value = state.scheme; onedSel.value = state.onedView; if (state.colorMax > 0) cmaxEl.value = state.colorMax;

    function show(sel, on) { const el = $(sel); if (el) el.style.display = on ? "" : "none"; }
    function updateControls() {
      const nG = getGOIs().length, nC = state.selectedCancers.length;
      const isMatrix = nG > 2 && nC === 1;
      const isOneD = nG === 2 && nC > 1;
      const isHeatmapColor = isMatrix || (isOneD && state.onedView === "heatmap");
      show("#cr-oned-wrap", isOneD);
      show("#cr-swap", isOneD && state.onedView === "heatmap");
      show("#cr-scheme-wrap", isHeatmapColor);
      show("#cr-cmax-wrap", isHeatmapColor);
    }

    // ---- 癌種 chips（複用 immuneCorr 互動）----
    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "cr-chip" + (state.selectedCancers.includes(code) ? " on" : "");
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
    $("#cr-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#cr-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#cr-def").addEventListener("click", () => { state.cancers = availCancers.slice(); commit(); renderCancers(); });
    $("#cr-call").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); });
    $("#cr-cnone").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); });

    methodSel.addEventListener("change", () => { state.corrMethod = methodSel.value; commit(); });
    schemeSel.addEventListener("change", () => { state.scheme = schemeSel.value; commit(); });
    onedSel.addEventListener("change", () => { state.onedView = onedSel.value; commit(); updateControls(); });
    cmaxEl.addEventListener("change", () => { const v = Number(cmaxEl.value); state.colorMax = (cmaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : Math.min(1, v); commit(); });
    $("#cr-swap").addEventListener("click", () => { state.swap = !state.swap; commit(); if (resultEl.querySelector("svg")) run(); });

    $("#cr-run").addEventListener("click", run);
    $("#cr-csv").addEventListener("click", exportCSV);
    $("#cr-svg").addEventListener("click", downloadSVG);
    $("#cr-png").addEventListener("click", downloadPNG);

    const methodName = () => state.corrMethod === "pearson" ? "Pearson" : "Spearman";

    // ---- 一對基因在某癌種的相關（pairwise-complete，丟 NaN）----
    function corrPair(gA, gB, cancer) {
      const aArr = geneVals.get(gA.rec.gene_id), bArr = geneVals.get(gB.rec.gene_id);
      const ps = byCancer[cancer] || [];
      const xs = [], ys = [];
      for (const p of ps) { const a = aArr[p.idx], b = bArr[p.idx]; if (isFinite(a) && isFinite(b)) { xs.push(a); ys.push(b); } }
      const n = xs.length;
      if (n < 4) return { n, state: "nodata", tip: `${gA.label} ~ ${gB.label} / ${cancer}: n=${n}` };
      const res = state.corrMethod === "pearson" ? pearsonr(xs, ys) : spearmanr(xs, ys);
      const ci = corrCI(res.r, n);
      const weak = n < 20;
      return { r: res.r, p: res.p, ciLow: ci.lo, ciHigh: ci.hi, n, xs, ys, state: weak ? "weak" : "ok", tip: `${gA.label} ~ ${gB.label} / ${cancer}: r=${res.r.toFixed(3)}, p=${res.p.toPrecision(2)}, n=${n}` };
    }

    async function run() {
      resultEl.innerHTML = ""; pairEl.innerHTML = ""; ["cr-csv", "cr-svg", "cr-png"].forEach(id => $("#" + id).style.display = "none");
      const cancers = state.cancers.filter(c => state.selectedCancers.includes(c));
      if (!cancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      const recs = [], unknown = [];
      getGOIs().forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) unknown.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (recs.length < 2) { statusEl.className = "status err"; statusEl.textContent = "Add at least 2 recognized genes." + (unknown.length ? ` Unrecognized: ${unknown.join(", ")}` : ""); return; }
      if (recs.length > 2 && cancers.length > 1) { statusEl.className = "status err"; statusEl.textContent = "Too many dimensions: with >2 genes pick a single cancer (for the matrix); with multiple cancers use exactly 2 genes."; return; }

      statusEl.className = "status"; statusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      try { geneVals = new Map(); await Promise.all(recs.map(async x => geneVals.set(x.rec.gene_id, await dataset.getGeneValues(x.rec)))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files (CORS?): " + e.message; return; }

      statusEl.textContent = "Computing…";
      const patients = patientsInScope(dataset, cancers);
      byCancer = {}; cancers.forEach(c => byCancer[c] = patients.filter(p => p.cancer === c));
      lastRecs = recs;

      let isScatter = false;
      if (recs.length === 2 && cancers.length === 1) { scatterWidget(recs[0], recs[1], cancers[0], resultEl); isScatter = true; }
      else if (recs.length === 2 && cancers.length > 1) { if (state.onedView === "bar") drawPairBar(recs[0], recs[1], cancers); else drawPairHeatmap(recs[0], recs[1], cancers); }
      else drawMatrix(recs, cancers[0]);

      const note = document.createElement("div"); note.className = "cr-legend";
      note.textContent = `${methodName()} · tumor only · per-cancer${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      // scatter 用 widget 自帶存檔；矩陣/1-D 才用主匯出鍵
      ["cr-csv", "cr-svg", "cr-png"].forEach(id => $("#" + id).style.display = isScatter ? "none" : "");
      statusEl.textContent = "Done.";
    }

    // ---- scatter widget（一對 × 一癌種；含控制區：X/Y 上下限、顏色、存圖/CSV）----
    function scatterWidget(gA, gB, cancer, hostEl) {
      const c = corrPair(gA, gB, cancer);
      if (c.state === "nodata") { hostEl.innerHTML = `<div class="status err">Too few paired samples (n=${c.n}).</div>`; return; }
      const xs = c.xs, ys = c.ys, baseName = `corr_${gA.label}_${gB.label}_${cancer}`;
      const fmt = v => Number(v.toPrecision(4));
      const dft = { xMin: fmt(Math.min(...xs)), xMax: fmt(Math.max(...xs)), yMin: fmt(Math.min(...ys)), yMax: fmt(Math.max(...ys)) };
      hostEl.innerHTML = `
        <div class="cr-pctl">
          <span class="cr-ptitle">${gA.label} vs ${gB.label} · ${cancer}</span>
          <span>X<input class="cr-pn" data-k="xMin" type="number" step="any"> – <input class="cr-pn" data-k="xMax" type="number" step="any"></span>
          <span>Y<input class="cr-pn" data-k="yMin" type="number" step="any"> – <input class="cr-pn" data-k="yMax" type="number" step="any"></span>
          <span>Opacity<input class="cr-popacity" type="range" min="0.1" max="1" step="0.05"></span>
          <span>Stats<select class="cr-plabel"><option value="tr">top-right</option><option value="tl">top-left</option></select></span>
          <button class="cr-mini" data-act="reset">Reset</button>
          <button class="cr-mini" data-act="svg">SVG</button>
          <button class="cr-mini" data-act="png">PNG</button>
          <button class="cr-mini" data-act="csv">CSV</button>
        </div>
        <div class="cr-pswatches"></div>
        <div class="cr-pplot"></div>`;
      const plotDiv = hostEl.querySelector(".cr-pplot"), swatchBox = hostEl.querySelector(".cr-pswatches");
      const opacityInp = hostEl.querySelector(".cr-popacity"), labelSel = hostEl.querySelector(".cr-plabel");
      opacityInp.value = lastScatterOpacity; labelSel.value = lastLabelPos;
      const setVals = o => Object.entries(o).forEach(([k, v]) => { hostEl.querySelector(`[data-k="${k}"]`).value = v; });
      setVals(dft);
      const numAt = k => { const v = Number(hostEl.querySelector(`[data-k="${k}"]`).value); return isFinite(v) ? v : null; };
      function renderSwatches() { swatchBox.innerHTML = PALETTE.map(col => `<span class="cr-sw${col === lastScatterColor ? " on" : ""}" data-col="${col}" style="background:${col}"></span>`).join(""); }
      function renderPlot() {
        plotDiv.innerHTML = corrScatterSVG(xs.map((x, i) => ({ x, y: ys[i] })), {
          xLabel: `${gA.label} (${unitOf(gA.rec)})`, yLabel: `${gB.label} (${unitOf(gB.rec)})`,
          r: c.r, pText: fmtP(c.p), n: c.n, title: `${gA.label} vs ${gB.label}`, subtitle: `${cancer} · ${methodName()}`,
          band: true, marginals: true, pointColor: lastScatterColor, pointOpacity: lastScatterOpacity, labelPos: lastLabelPos,
          xMin: numAt("xMin"), xMax: numAt("xMax"), yMin: numAt("yMin"), yMax: numAt("yMax"),
        });
      }
      renderSwatches(); renderPlot();
      hostEl.querySelectorAll(".cr-pn").forEach(el => el.addEventListener("input", renderPlot));
      swatchBox.addEventListener("click", e => { const sw = e.target.closest(".cr-sw"); if (!sw) return; lastScatterColor = sw.dataset.col; renderSwatches(); renderPlot(); });
      opacityInp.addEventListener("input", () => { lastScatterOpacity = Number(opacityInp.value); renderPlot(); });
      labelSel.addEventListener("change", () => { lastLabelPos = labelSel.value; renderPlot(); });
      hostEl.querySelector('[data-act="reset"]').addEventListener("click", () => { setVals(dft); renderPlot(); });
      hostEl.querySelector('[data-act="svg"]').addEventListener("click", () => dlSVGEl(plotDiv.querySelector("svg"), baseName));
      hostEl.querySelector('[data-act="png"]').addEventListener("click", () => dlPNGEl(plotDiv.querySelector("svg"), baseName));
      hostEl.querySelector('[data-act="csv"]').addEventListener("click", () => { let csv = `${csvq(gA.label)},${csvq(gB.label)}\n`; for (let i = 0; i < xs.length; i++) csv += `${xs[i].toFixed(4)},${ys[i].toFixed(4)}\n`; dl(new Blob([csv], { type: "text/csv" }), baseName + ".csv"); });
    }

    // ---- >2 基因 × 1 癌種：上三角矩陣（點格子 → scatter）----
    function drawMatrix(recs, cancer) {
      lastSVGName = `corr_matrix_${cancer}`;
      const N = recs.length;
      const cells = recs.map(() => new Array(N).fill(null));
      const pflat = [];
      for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
        const c = corrPair(recs[i], recs[j], cancer);
        if (c.state === "nodata") cells[i][j] = { state: "nodata", tip: c.tip };
        else { cells[i][j] = { value: c.r, state: c.state, stars: "", tip: c.tip }; pflat.push({ i, j, p: c.p }); }
      }
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = cells[x.i][x.j]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      lastResult = { type: "matrix", labels: recs.map(r => r.label), cells };
      const colorMax = state.colorMax > 0 ? state.colorMax : 1;
      resultEl.innerHTML = corrMatrixSVG(recs.map(r => r.label), cells, { colorMax, scheme: state.scheme, legendLabel: `${methodName()} r`, caption: `Gene correlation · ${cancer}` });
      // 點上三角格子 → 該對 scatter（顯示在下方 pair 區）
      const svg = resultEl.querySelector("svg");
      if (svg) svg.addEventListener("click", e => {
        const t = e.target.closest && e.target.closest("[data-i]"); if (!t) return;
        scatterWidget(recs[+t.dataset.i], recs[+t.dataset.j], cancer, pairEl);
        pairEl.scrollIntoView({ behavior: "smooth", block: "nearest" });
      });
    }

    // ---- 2 基因 × 多癌種：r across cancers ----
    function pairAcrossCancers(gA, gB, cancers) {
      return cancers.map(c => ({ cancer: c, c: corrPair(gA, gB, c) }));
    }
    function drawPairBar(gA, gB, cancers) {
      lastSVGName = `corr_bar_${gA.label}_${gB.label}`;
      const rows = pairAcrossCancers(gA, gB, cancers);
      const pv = rows.filter(x => x.c.state !== "nodata").map(x => x.c.p);
      const q = benjaminiHochberg(pv); let qi = 0;
      rows.forEach(x => { if (x.c.state !== "nodata") { const qq = q[qi++]; x.c.q = qq; const st = pStars(qq); x.c.stars = st === "ns" ? "" : st; x.c.tip += `, q=${qq.toPrecision(2)}`; } });
      lastResult = { type: "oned", gA, gB, rows };
      const items = rows.map(x => x.c.state === "nodata" ? { label: x.cancer, state: "nodata" } : { label: x.cancer, r: x.c.r, ciLow: x.c.ciLow, ciHigh: x.c.ciHigh, stars: x.c.stars, state: x.c.state, tip: x.c.tip });
      resultEl.innerHTML = corrBarSVG(items, { caption: `${gA.label} vs ${gB.label}`, yLabel: `${methodName()} r`, scheme: state.scheme });
    }
    function drawPairHeatmap(gA, gB, cancers) {
      lastSVGName = `corr_1d_${gA.label}_${gB.label}`;
      const rows = pairAcrossCancers(gA, gB, cancers);
      const pflat = []; rows.forEach((x, r) => { if (x.c.state !== "nodata") pflat.push({ r, p: x.c.p }); });
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const c = rows[x.r].c; c.q = q[k]; const st = pStars(q[k]); c.stars = st === "ns" ? "" : st; c.tip += `, q=${q[k].toPrecision(2)}`; });
      lastResult = { type: "oned", gA, gB, rows };
      const cellObjs = rows.map(x => x.c.state === "nodata" ? { state: "nodata", tip: x.c.tip } : { value: x.c.r, state: x.c.state, stars: x.c.stars, tip: x.c.tip });
      const colorMax = state.colorMax > 0 ? state.colorMax : 1;
      const caption = `${gA.label} vs ${gB.label}\n${methodName()} r across cancers`;   // A vs B 換行、不用 ·
      const cancerLabels = rows.map(x => x.cancer);
      // 標題已標 GOIs → 不重複欄/列標籤；swap 切直式(癌種為列) / 橫式(癌種為欄)
      if (state.swap) resultEl.innerHTML = heatmapSVG([""], cancerLabels, [cellObjs], { scheme: state.scheme, colorMax, legendLabel: `${methodName()} r`, caption });
      else resultEl.innerHTML = heatmapSVG(cancerLabels, [""], cellObjs.map(c => [c]), { scheme: state.scheme, colorMax, legendLabel: `${methodName()} r`, caption });
    }

    // ---- 匯出 CSV（缺值/不可估計 → NA）----
    function exportCSV() {
      if (!lastResult) return;
      let csv = "";
      if (lastResult.type === "matrix") {
        const { labels, cells } = lastResult, N = labels.length;       // 對稱方陣：對角=1、上三角填、下三角鏡射、nodata→NA
        csv = "," + labels.map(csvq).join(",") + "\n";
        for (let i = 0; i < N; i++) {
          const row = [csvq(labels[i])];
          for (let j = 0; j < N; j++) {
            if (i === j) { row.push("1"); continue; }
            const c = i < j ? cells[i][j] : cells[j][i];
            row.push(c && c.state !== "nodata" && isFinite(c.value) ? c.value.toFixed(4) : "NA");
          }
          csv += row.join(",") + "\n";
        }
      } else if (lastResult.type === "oned") {
        csv = "cancer,r,p,q,n\n";
        lastResult.rows.forEach(x => { const c = x.c; csv += c.state === "nodata" ? `${csvq(x.cancer)},NA,NA,NA,${c.n || 0}\n` : `${csvq(x.cancer)},${c.r.toFixed(4)},${c.p.toExponential(4)},${c.q != null ? c.q.toExponential(4) : "NA"},${c.n}\n`; });
      } else if (lastResult.type === "scatter") {
        const { gA, gB, xs, ys } = lastResult;
        csv = `${csvq(gA.label)},${csvq(gB.label)}\n`;
        for (let i = 0; i < xs.length; i++) csv += `${xs[i].toFixed(4)},${ys[i].toFixed(4)}\n`;
      } else return;
      dl(new Blob([csv], { type: "text/csv" }), lastSVGName + ".csv");
    }

    // ---- 下載（主結果區）----
    function firstSVG() { return resultEl.querySelector("svg") || pairEl.querySelector("svg"); }
    function dl(blob, name) { const u = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = u; a.download = name; a.click(); URL.revokeObjectURL(u); }
    function dlSVGEl(svg, name) { if (!svg) return; dl(new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml;charset=utf-8" }), name + ".svg"); }
    function dlPNGEl(svg, name) {
      if (!svg) return;
      const vb = svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal : null;
      const W = vb && vb.width ? vb.width : svg.clientWidth, H = vb && vb.height ? vb.height : svg.clientHeight;
      const scale = Math.max(2, Math.ceil(2400 / W));   // 不論畫面縮多小，輸出固定高解析度
      const xml = new XMLSerializer().serializeToString(svg);
      const img = new Image();
      img.onload = () => { const c = document.createElement("canvas"); c.width = Math.round(W * scale); c.height = Math.round(H * scale); const ctx = c.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, c.width, c.height); ctx.drawImage(img, 0, 0, c.width, c.height); c.toBlob(b => dl(b, name + ".png"), "image/png"); };
      img.onerror = () => alert("PNG export failed (try SVG).");
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(xml);
    }
    function downloadSVG() { dlSVGEl(firstSVG(), lastSVGName); }
    function downloadPNG() { dlPNGEl(firstSVG(), lastSVGName); }

    renderCancers();
  },
};
