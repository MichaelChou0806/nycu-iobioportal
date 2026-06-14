// =====================================================================
// analyses/survival.js —「Survival (KM)」（OS）
// ---------------------------------------------------------------------
// 只取 Tumor。GOIs 與 clinicalOverview 共用（core/gois.js）。三種切法
//   median / tertile / quartile（後兩者丟中間組）。隨訪月數可截斷。
//   - 單基因 + 一個癌種 → KM 曲線 + log-rank p + Cox HR(95%CI) + 每組 n
//   - 單基因 + 多癌種：可選
//        · Per-cancer → HR heatmap（紅/綠或紅/藍，★=log-rank FDR）
//        · Pooled     → 一張 KM（每癌種內各自切 High/Low，按癌種分層的 log-rank/Cox）
//   - 多基因 → HR heatmap
// 可調色階上下限、紅綠/紅藍配色（色盲友善）。KM 可匯出 Prism 用 CSV。
// 統計核心見 stats.js（已驗證）。狀態自存（癌種/切法等，GOIs 走共享）。
// =====================================================================

import { patientsInScope } from "../core/dimensions.js";
import { kaplanMeier, logRank, coxPH1, logRankStratified, coxPH1Stratified, benjaminiHochberg, pStars } from "../core/stats.js";
import { kmCurveSVG, heatmapSVG } from "../core/plots.js";
import { getGOIs, setGOIs, onGOIsChanged, parseGenes } from "../core/gois.js";

const DAYS_PER_MONTH = 30.4375;
const HIGH_COLOR = "#ef4444", LOW_COLOR = "#3b82f6";
const SHORT = { median: "median", tertile: "tertile", quartile: "quartile" };

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .sv-sec{margin-bottom:18px}
    .sv-sum{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;font-weight:600;cursor:pointer;padding:4px 0}
    details.sv-sec>summary{margin-bottom:6px}
    .sv-h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .sv-genes{width:100%;min-height:56px;padding:9px 11px;border:1px solid var(--line);border-radius:7px;font:13px/1.5 monospace;resize:vertical}
    .sv-toolbar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;align-items:center}
    .sv-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .sv-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .sv-chips{display:flex;flex-wrap:wrap;gap:6px}
    .sv-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .sv-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .sv-chip.dragging{opacity:.4}
    .sv-note{font-size:11.5px;color:var(--muted)}
    .sv-result{overflow-x:auto}
    .sv-result svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .sv-legend{font-size:11.5px;color:var(--muted);margin-top:6px}
    .sv-ctl{font-size:11.5px;color:var(--muted)}
    .sv-ctl select,.sv-ctl input{padding:5px 8px;font-size:12px;margin-left:4px}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

const LS = "tcga-tool:survival";
function loadState() { try { return JSON.parse(localStorage.getItem(LS)) || {}; } catch (e) { return {}; } }
function saveState(st) { try { localStorage.setItem(LS, JSON.stringify(st)); } catch (e) {} }
function moveInArray(arr, fromId, toId) {
  if (fromId == null || fromId === toId) return arr.slice();
  const a = arr.slice(); const fi = a.indexOf(fromId); if (fi < 0) return a;
  a.splice(fi, 1); let ti = a.indexOf(toId); if (ti < 0) ti = a.length; a.splice(ti, 0, fromId); return a;
}
function quantileSorted(sorted, f) { const p = f * (sorted.length - 1), lo = Math.floor(p), hi = Math.ceil(p); return sorted[lo] + (sorted[hi] - sorted[lo]) * (p - lo); }
// 切法分組：1=High, 0=Low, -1=丟棄
function splitGroups(vals, method) {
  const sorted = [...vals].sort((a, b) => a - b);
  const frac = method === "tertile" ? 1 / 3 : method === "quartile" ? 1 / 4 : 0.5;
  const lowCut = quantileSorted(sorted, frac), highCut = quantileSorted(sorted, 1 - frac);
  return vals.map(v => method === "median" ? (v > highCut ? 1 : 0) : (v <= lowCut ? 0 : (v > highCut ? 1 : -1)));
}
function splitLabel(m) { return m === "tertile" ? "tertile (top/bottom 1/3)" : m === "quartile" ? "quartile (Q1 vs Q4)" : "median split"; }
function fmtP(p) { return p < 0.001 ? "< 0.001" : p.toFixed(3); }

export const survival = {
  id: "survival",
  name: "Survival (KM)",

  async mount(container, { dataset }) {
    injectStyles();
    const availCancers = dataset.cancers.map(c => c.code);
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));
    const saved = loadState();
    let state = {
      cancers: (Array.isArray(saved.cancers) ? saved.cancers.filter(c => availCancers.includes(c)) : []).concat(availCancers.filter(c => !(saved.cancers || []).includes(c))),
      selectedCancers: Array.isArray(saved.selectedCancers) ? saved.selectedCancers.filter(c => availCancers.includes(c)) : [],
      split: saved.split || "median",
      months: saved.months != null ? saved.months : 0,
      view: saved.view || "percancer",     // 單基因多癌種：percancer | pooled
      scheme: saved.scheme || "rg",          // rg=紅綠, rb=紅藍(色盲友善)
      colorMax: saved.colorMax != null ? saved.colorMax : 0,  // 0=auto
    };
    let dragCancer = null, lastSVGName = "km", lastKM = null;  // lastKM 供 Prism 匯出

    container.innerHTML = `
      <div class="card">
        <div class="sv-sec">
          <h3 class="sv-h3">Genes (shared with Clinical Overview)</h3>
          <textarea id="sv-genes" class="sv-genes"></textarea>
        </div>

        <details class="sv-sec" open>
          <summary class="sv-sum">Cancers — drag to reorder · pick ONE for a KM curve, many for heatmap/pooled</summary>
          <div class="sv-toolbar" style="margin-top:8px">
            <button class="sv-mini" id="sv-alpha">Sort A–Z</button>
            <button class="sv-mini" id="sv-byn">Sort by N</button>
            <button class="sv-mini" id="sv-def">Default order</button>
            <span class="sv-sep"></span>
            <button class="sv-mini" id="sv-all">Select all</button>
            <button class="sv-mini" id="sv-none">Select none</button>
          </div>
          <div id="sv-cancers" class="sv-chips"></div>
        </details>

        <div class="sv-sec">
          <div class="sv-toolbar">
            <label class="sv-ctl">Split
              <select id="sv-split">
                <option value="median">Median (50/50)</option>
                <option value="tertile">Tertile (top/bottom 1/3)</option>
                <option value="quartile">Quartile (Q4 vs Q1)</option>
              </select></label>
            <label class="sv-ctl">Months<input type="number" id="sv-months" min="0" step="1" style="width:72px" placeholder="all"></label>
            <label class="sv-ctl" id="sv-view-wrap">Multi-cancer
              <select id="sv-view">
                <option value="percancer">Per-cancer (HR heatmap)</option>
                <option value="pooled">Pooled (one KM, stratified)</option>
              </select></label>
            <label class="sv-ctl" id="sv-scheme-wrap">Colors
              <select id="sv-scheme">
                <option value="rg">Red–Green</option>
                <option value="rb">Red–Blue (CB-friendly)</option>
              </select></label>
            <label class="sv-ctl" id="sv-cmax-wrap">Scale max<input type="number" id="sv-cmax" min="0" step="0.1" style="width:64px" placeholder="auto"></label>
            <button id="sv-run">Run</button>
            <button class="sv-mini" id="sv-prism" style="display:none">Export Prism CSV</button>
            <button class="sv-mini" id="sv-svg" style="display:none">Download SVG</button>
            <button class="sv-mini" id="sv-png" style="display:none">Download PNG</button>
          </div>
          <div id="sv-status" class="status"></div>
          <div id="sv-result" class="sv-result"></div>
        </div>
        <div class="sv-note">Tumor only · endpoint OS. Tertile/Quartile drop the middle group (per-group n shown). Pooled splits High/Low within each cancer and uses cancer-stratified log-rank/Cox.</div>
      </div>`;

    const $ = s => container.querySelector(s);
    const chipBox = $("#sv-cancers"), genesEl = $("#sv-genes"), statusEl = $("#sv-status"), resultEl = $("#sv-result");
    const splitSel = $("#sv-split"), monthsEl = $("#sv-months"), viewSel = $("#sv-view"), schemeSel = $("#sv-scheme"), cmaxEl = $("#sv-cmax");

    function commit() { saveState(state); }
    // GOIs 共享
    genesEl.value = getGOIs().join(", ");
    genesEl.addEventListener("input", () => setGOIs(parseGenes(genesEl.value)));
    const off = onGOIsChanged(list => { if (document.activeElement !== genesEl) genesEl.value = list.join(", "); });
    // 控制項初值
    splitSel.value = state.split; if (state.months > 0) monthsEl.value = state.months;
    viewSel.value = state.view; schemeSel.value = state.scheme; if (state.colorMax > 0) cmaxEl.value = state.colorMax;

    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "sv-chip" + (state.selectedCancers.includes(code) ? " on" : "");
        chip.dataset.code = code; chip.textContent = code; chip.draggable = true;
        chip.addEventListener("click", () => {
          if (state.selectedCancers.includes(code)) { state.selectedCancers = state.selectedCancers.filter(c => c !== code); chip.classList.remove("on"); }
          else { state.selectedCancers.push(code); chip.classList.add("on"); }
          commit();
        });
        chip.addEventListener("dragstart", () => { dragCancer = code; chip.classList.add("dragging"); });
        chip.addEventListener("dragend", () => chip.classList.remove("dragging"));
        chip.addEventListener("dragover", e => e.preventDefault());
        chip.addEventListener("drop", e => { e.preventDefault(); state.cancers = moveInArray(state.cancers, dragCancer, code); commit(); renderCancers(); });
        chipBox.appendChild(chip);
      });
    }
    $("#sv-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#sv-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#sv-def").addEventListener("click", () => { state.cancers = availCancers.slice(); commit(); renderCancers(); });
    $("#sv-all").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); });
    $("#sv-none").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); });
    splitSel.addEventListener("change", () => { state.split = splitSel.value; commit(); });
    monthsEl.addEventListener("change", () => { const v = Number(monthsEl.value); state.months = (monthsEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); });
    viewSel.addEventListener("change", () => { state.view = viewSel.value; commit(); });
    schemeSel.addEventListener("change", () => { state.scheme = schemeSel.value; commit(); });
    cmaxEl.addEventListener("change", () => { const v = Number(cmaxEl.value); state.colorMax = (cmaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); });

    // 取存活資料（只取 tumor、套月數截斷）：回傳 {expr,e,tm}[]
    function survivalRows(vals, patients) {
      const cut = state.months > 0 ? state.months : Infinity; const out = [];
      for (const p of patients) {
        const expr = vals[p.idx], os = Number(p.clin["OS"]), t = Number(p.clin["OS.time"]);
        if (!isFinite(expr) || !isFinite(os) || !isFinite(t) || t < 0) continue;
        let tm = t / DAYS_PER_MONTH, e = os ? 1 : 0;
        if (tm > cut) { e = 0; tm = cut; }
        out.push({ expr, e, tm });
      }
      return out;
    }
    function highLow(rows) {
      const g = splitGroups(rows.map(r => r.expr), state.split); const high = [], low = [];
      rows.forEach((r, i) => { if (g[i] === 1) high.push(r); else if (g[i] === 0) low.push(r); });
      return { high, low };
    }

    // ---- RUN ----
    $("#sv-run").addEventListener("click", run);
    $("#sv-prism").addEventListener("click", exportPrism);
    $("#sv-svg").addEventListener("click", downloadSVG);
    $("#sv-png").addEventListener("click", downloadPNG);

    async function run() {
      resultEl.innerHTML = ""; lastKM = null;
      ["sv-prism", "sv-svg", "sv-png"].forEach(id => $("#" + id).style.display = "none");
      const cancers = state.cancers.filter(c => state.selectedCancers.includes(c));
      if (!cancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      const recs = []; const unknown = [];
      getGOIs().forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) unknown.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (!recs.length) { statusEl.className = "status err"; statusEl.textContent = "No recognized genes." + (unknown.length ? ` Unrecognized: ${unknown.join(", ")}` : ""); return; }

      statusEl.className = "status"; statusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      let geneVals;
      try { geneVals = new Map(); await Promise.all(recs.map(async x => geneVals.set(x.rec.gene_id, await dataset.getGeneValues(x.rec)))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files (CORS?): " + e.message; return; }

      statusEl.textContent = "Computing…";
      const patients = patientsInScope(dataset, cancers);
      const single = recs.length === 1;
      let isKM = false;
      if (single && cancers.length === 1) { drawKM(recs[0], cancers[0], geneVals, patients); isKM = true; }
      else if (single && state.view === "pooled") { drawPooledKM(recs[0], cancers, geneVals, patients); isKM = true; }
      else drawHeatmap(recs, cancers, geneVals, patients);

      const note = document.createElement("div"); note.className = "sv-legend";
      note.textContent = isKM
        ? `Endpoint OS · ${splitLabel(state.split)}${state.months > 0 ? ` · capped ${state.months} mo` : ""}${lastKM && lastKM.stratified ? " · cancer-stratified" : ""}${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`
        : `★ log-rank FDR q<0.05 ★★<0.01 ★★★<0.001 ★★★★<0.0001 | red=HR>1 (worse) ${state.scheme === "rb" ? "blue" : "green"}=HR<1 (better) | grey=no data, faded=too few${unknown.length ? ` | Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      if (isKM) $("#sv-prism").style.display = "";
      $("#sv-svg").style.display = ""; $("#sv-png").style.display = "";
      statusEl.textContent = "Done.";
    }

    // 把 High/Low 兩組畫成 KM；stratLabel 有值時用分層統計
    function renderKM(geneLabel, scopeText, high, low, strata, captionScope) {
      if (!high.length || !low.length) { resultEl.innerHTML = `<div class="status err">Not enough samples to split (High n=${high.length}, Low n=${low.length}).</div>`; return false; }
      const kmHigh = kaplanMeier(high.map(r => r.tm), high.map(r => r.e));
      const kmLow = kaplanMeier(low.map(r => r.tm), low.map(r => r.e));
      const allTm = high.concat(low).map(r => r.tm), allE = high.concat(low).map(r => r.e), allG = high.map(() => 1).concat(low.map(() => 0));
      let lr, cox, stratified = false;
      if (strata) { const st = strata; lr = logRankStratified(allTm, allE, allG, st); cox = coxPH1Stratified(allTm, allE, allG, st); stratified = true; }
      else { lr = logRank(allTm, allE, allG); cox = coxPH1(allTm, allE, allG); }
      const curves = [{ label: "High", color: HIGH_COLOR, km: kmHigh, n: high.length }, { label: "Low", color: LOW_COLOR, km: kmLow, n: low.length }];
      const hrText = isFinite(cox.hr) ? `HR (High vs Low${stratified ? ", stratified" : ""}) = ${cox.hr.toFixed(2)} (95% CI ${cox.ciLow.toFixed(2)}–${cox.ciHigh.toFixed(2)})` : "HR = n/a";
      resultEl.innerHTML = kmCurveSVG(curves, {
        caption: captionScope,
        xLabel: state.months > 0 ? `Months (capped ${state.months})` : "Months",
        pText: `Log-rank p = ${fmtP(lr.p)}${stratified ? " (stratified)" : ""}`, hrText,
        xMax: state.months > 0 ? state.months : undefined,
      });
      lastKM = { geneLabel, high, low, stratified, scope: scopeText };
      return true;
    }

    // 單基因單癌種
    function drawKM(geneRec, cancer, geneVals, patients) {
      lastSVGName = `KM_${geneRec.label}_${cancer}`;
      const vals = geneVals.get(geneRec.rec.gene_id);
      const { high, low } = highLow(survivalRows(vals, patients.filter(p => p.cancer === cancer)));
      renderKM(geneRec.label, `${geneRec.label} ${cancer}`, high, low, null, `${geneRec.label} — ${cancer} · OS · ${SHORT[state.split]}`);
    }

    // 單基因多癌種 pooled：每癌種內各自切，合併，按癌種分層
    function drawPooledKM(geneRec, cancers, geneVals, patients) {
      lastSVGName = `KM_${geneRec.label}_pooled`;
      const vals = geneVals.get(geneRec.rec.gene_id);
      // 逐癌種內各自切 High/Low，並同步記錄每筆的分層（癌種）
      const high = [], low = [], hStrata = [], lStrata = [];
      cancers.forEach(c => {
        const { high: h, low: l } = highLow(survivalRows(vals, patients.filter(p => p.cancer === c)));
        h.forEach(r => { high.push(r); hStrata.push(c); });
        l.forEach(r => { low.push(r); lStrata.push(c); });
      });
      const strataAll = hStrata.concat(lStrata);  // 對齊 high.concat(low)
      renderKM(geneRec.label, `${geneRec.label} pooled`, high, low, strataAll,
        `${geneRec.label} · pooled ${cancers.length} cancers (stratified) · OS · ${SHORT[state.split]}`);
    }

    // 多基因/多癌種：HR heatmap
    function drawHeatmap(recs, cancers, geneVals, patients) {
      lastSVGName = "HR_heatmap";
      const rows = recs.map(r => r.label), cols = cancers;
      const byCancer = {}; cancers.forEach(c => byCancer[c] = patients.filter(p => p.cancer === c));
      const cells = []; const pflat = [];
      recs.forEach((r, ri) => {
        const vals = geneVals.get(r.rec.gene_id); const rowCells = [];
        cancers.forEach((c, ci) => {
          const { high, low } = highLow(survivalRows(vals, byCancer[c]));
          const evHigh = high.reduce((s, x) => s + x.e, 0), evLow = low.reduce((s, x) => s + x.e, 0);
          let cell;
          if (!high.length || !low.length) cell = { state: "nodata", tip: `${r.label} / ${c}: cannot split (n=${high.length}/${low.length})` };
          else if (evHigh + evLow === 0) cell = { state: "nodata", tip: `${r.label} / ${c}: no events` };
          else {
            const allTm = high.concat(low).map(x => x.tm), allE = high.concat(low).map(x => x.e), allG = high.map(() => 1).concat(low.map(() => 0));
            const lr = logRank(allTm, allE, allG), cox = coxPH1(allTm, allE, allG);
            const log2hr = isFinite(cox.hr) && cox.hr > 0 ? Math.log2(cox.hr) : 0;
            const weak = high.length < 10 || low.length < 10 || (evHigh + evLow) < 10;
            cell = { value: log2hr, state: weak ? "weak" : "ok", stars: "", tip: `${r.label} / ${c}: HR=${cox.hr.toFixed(2)} (${cox.ciLow.toFixed(2)}–${cox.ciHigh.toFixed(2)}), n=${high.length}/${low.length}, events=${evHigh + evLow}, p=${lr.p.toPrecision(2)}` };
            if (!weak) pflat.push({ ri, ci, p: lr.p });
          }
          rowCells.push(cell);
        });
        cells.push(rowCells);
      });
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = cells[x.ri][x.ci]; cell.q = q[k]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      let maxAbs = 0; cells.forEach(row => row.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) maxAbs = Math.max(maxAbs, Math.abs(c.value)); }));
      const colorMax = state.colorMax > 0 ? state.colorMax : Math.min(4, Math.max(0.5, Math.ceil(maxAbs * 10) / 10));
      resultEl.innerHTML = heatmapSVG(rows, cols, cells, {
        colorMax, scheme: state.scheme, legendLabel: "log2 HR",
        caption: `Hazard ratio (High vs Low) · OS · ${SHORT[state.split]}${state.months > 0 ? ` · ${state.months} mo` : ""}`,
      });
    }

    // ---- 匯出：Prism 生存格式（Month, GOI Low, GOI High；對不到該組的留空）----
    function exportPrism() {
      if (!lastKM) return;
      const g = lastKM.geneLabel;
      const rows = [];
      lastKM.low.forEach(r => rows.push({ m: r.tm, low: r.e, high: "" }));
      lastKM.high.forEach(r => rows.push({ m: r.tm, low: "", high: r.e }));
      rows.sort((a, b) => a.m - b.m);
      let csv = `Month,${g} Low,${g} High\n` + rows.map(r => `${r.m.toFixed(3)},${r.low},${r.high}`).join("\n") + "\n";
      dl(new Blob([csv], { type: "text/csv" }), `prism_KM_${g}_${lastKM.scope.replace(/\s+/g, "_")}.csv`);
    }
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
