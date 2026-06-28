// =====================================================================
// analyses/coxRegression.js —「Cox Regression」（univariate；OS）
// ---------------------------------------------------------------------
// 對單一癌種，逐因子做 univariate Cox（每因子一個 HR），畫成 forest plot。
//   - 基因（與其他分析共用 GOIs）：High vs Low（median/tertile/quartile）
//   - 臨床因子（dimensions.tcga.json）：advanced vs baseline（無資料自動灰掉）
//   - vital 排除：它就是 OS 事件本身（用它預測死亡 = 完全分離）
// FDR 星號跨因子校正。HR=advanced/High 相對 baseline/Low（沿用 A=1,B=0 編碼）。
// 統計核心見 stats.js（coxPH1 已驗證）。狀態自存；GOIs 走共享。
// 之後會在同檔加 multivariate（用 stats.js 的 coxPH 多共變量引擎）。
// =====================================================================

import { patientsInScope, loadDimensions, analyzeDimension, classify } from "../core/dimensions.js";
import { coxPH1, coxPH, benjaminiHochberg, pStars } from "../core/stats.js";
import { forestSVG, heatmapSVG } from "../core/plots.js";
import { getGOIs, setGOIs, onGOIsChanged, parseGenes } from "../core/gois.js";

const DAYS_PER_MONTH = 30.4375;

// 存活 endpoint。實際顯示哪些由 clinical 有無對應欄位決定（資料驅動，OSCC 自動沿用）。
const ENDPOINTS = [
  { id: "OS", event: "OS", time: "OS.time", label: "OS" },
  { id: "DSS", event: "DSS", time: "DSS.time", label: "DSS" },
  { id: "DFI", event: "DFI", time: "DFI.time", label: "DFI" },
  { id: "PFI", event: "PFI", time: "PFI.time", label: "PFI" },
];

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .cx-sec{margin-bottom:18px}
    .cx-h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .cx-genes{width:100%;min-height:52px;padding:9px 11px;border:1px solid var(--line);border-radius:7px;font:13px/1.5 monospace;resize:vertical}
    .cx-toolbar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;align-items:center}
    .cx-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .cx-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .cx-ctl{font-size:11.5px;color:var(--muted)}
    .cx-ctl select,.cx-ctl input{padding:5px 8px;font-size:12px;margin-left:4px}
    .cx-seg{display:inline-flex;border:1px solid var(--line);border-radius:7px;overflow:hidden;margin-left:6px}
    .cx-seg-btn{padding:6px 14px;font-size:12.5px;font-weight:500;background:#fff;color:var(--ink);border:none;border-right:1px solid var(--line);border-radius:0;cursor:pointer}
    .cx-seg-btn:last-child{border-right:none}
    .cx-seg-btn.active{background:var(--accent);color:#fff}
    .cx-chips{display:flex;flex-wrap:wrap;gap:6px}
    .cx-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .cx-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .cx-chip.dragging{opacity:.4}
    .cx-dims{display:grid;grid-template-columns:repeat(auto-fill,minmax(210px,1fr));gap:8px}
    .cx-dim{display:flex;align-items:center;gap:8px;padding:7px 11px;border:1px solid var(--line);border-radius:6px;font-size:12.5px;cursor:pointer;user-select:none;min-width:0}
    .cx-dim.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .cx-dim.off{opacity:.45;cursor:not-allowed}
    .cx-dim.warn{border-color:#f0c089}
    .cx-dim.on small{color:#dbeafe}
    .cx-dag{color:#b45309;font-weight:700;font-size:11px;margin-left:2px}
    .cx-dim input{appearance:auto;width:15px;height:15px;min-width:0;padding:0;border:0;margin:0;cursor:pointer;flex:none}
    .cx-dim-name{flex:none;font-weight:500}
    .cx-dim small{flex:1 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;opacity:.75;font-size:10.5px}
    .cx-note{font-size:11.5px;color:var(--muted);line-height:1.5}
    .cx-result{overflow-x:auto}
    .cx-result svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .cx-legend{font-size:11.5px;color:var(--muted);margin-top:6px;line-height:1.5}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

const LS = "tcga-tool:cox";
function loadState() { try { return JSON.parse(localStorage.getItem(LS)) || {}; } catch (e) { return {}; } }
function saveState(st) { try { localStorage.setItem(LS, JSON.stringify(st)); } catch (e) {} }

function quantileSorted(sorted, f) { const p = f * (sorted.length - 1), lo = Math.floor(p), hi = Math.ceil(p); return sorted[lo] + (sorted[hi] - sorted[lo]) * (p - lo); }
// 切法分組：1=High, 0=Low, -1=丟棄（與 survival.js 同邏輯）
function splitGroups(vals, method) {
  const sorted = [...vals].sort((a, b) => a - b);
  const frac = method === "tertile" ? 1 / 3 : method === "quartile" ? 1 / 4 : 0.5;
  const lowCut = quantileSorted(sorted, frac), highCut = quantileSorted(sorted, 1 - frac);
  return vals.map(v => method === "median" ? (v > highCut ? 1 : 0) : (v <= lowCut ? 0 : (v > highCut ? 1 : -1)));
}
// 維度的 baseline / advanced 顯示標籤
function labelsOf(dim) {
  if (dim.type === "binary") {
    if (dim.numericSplit) return [dim.numericSplit.baselineLabel, dim.numericSplit.advancedLabel];
    return [dim.baseline.label, dim.advanced.label];
  }
  const base = dim.levels.filter(l => l.default === "baseline").map(l => l.id);
  const adv = dim.levels.filter(l => l.default === "advanced").map(l => l.id);
  return [base.join("/"), adv.join("/")];
}
function fmtP(p) { return p < 1e-4 ? "p<0.0001" : "p=" + (p < 0.001 ? p.toExponential(1) : p.toFixed(3)); }
function moveInArray(arr, fromId, toId) {
  if (fromId == null || fromId === toId) return arr.slice();
  const a = arr.slice(); const fi = a.indexOf(fromId); if (fi < 0) return a;
  a.splice(fi, 1); let ti = a.indexOf(toId); if (ti < 0) ti = a.length; a.splice(ti, 0, fromId); return a;
}

export const coxRegression = {
  id: "coxRegression",
  name: "Cox Regression",

  async mount(container, { dataset }) {
    injectStyles();
    const availCancers = dataset.cancers.map(c => c.code);
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));

    const dimUrl = (dataset.config && dataset.config.dimensionsUrl) || "config/dimensions.tcga.json";
    let DEF;
    try { DEF = await loadDimensions(dimUrl); }
    catch (e) { container.innerHTML = `<div class="card"><div class="status err">Failed to load clinical dimensions: ${e.message}</div></div>`; return; }
    // vital 排除（= 存活事件本身，當共變量會完全分離）。
    // recurrence 不硬排除：post-baseline，預設不選但可手動勾並標示。
    const DISCOURAGED = { recurrence: "post-baseline event" };
    const dims = (DEF.dimensions || []).filter(d => d.id !== "vital");

    // 可用 endpoint：clinical 同時有 event 與 time 欄位才列入（資料驅動）
    const clinFields = new Set(dataset.clinicalFields || (dataset.clinical.size ? Object.keys(dataset.clinical.values().next().value) : []));
    const endpoints = ENDPOINTS.filter(ep => clinFields.has(ep.event) && clinFields.has(ep.time));
    if (!endpoints.length) endpoints.push(ENDPOINTS[0]);   // 退路：至少有 OS

    const saved = loadState();
    let state = {
      endpoint: endpoints.some(ep => ep.id === saved.endpoint) ? saved.endpoint : endpoints[0].id,
      cancers: (Array.isArray(saved.cancers) ? saved.cancers.filter(c => availCancers.includes(c)) : []).concat(availCancers.filter(c => !(saved.cancers || []).includes(c))),
      selectedCancers: Array.isArray(saved.selectedCancers) ? saved.selectedCancers.filter(c => availCancers.includes(c)) : ((saved.cancer && availCancers.includes(saved.cancer)) ? [saved.cancer] : (availCancers.includes("HNSC") ? ["HNSC"] : availCancers.slice(0, 1))),
      selectedDims: Array.isArray(saved.selectedDims) ? saved.selectedDims.filter(id => dims.some(d => d.id === id)) : dims.filter(d => !DISCOURAGED[d.id]).map(d => d.id),
      split: ["median", "tertile", "quartile"].includes(saved.split) ? saved.split : "median",
      months: saved.months != null ? saved.months : 0,
      model: saved.model === "multi" ? "multi" : "uni",
      swap: !!saved.swap,
      scheme: saved.scheme === "rb" ? "rb" : "rg",
    };
    const currentEndpoint = () => endpoints.find(ep => ep.id === state.endpoint) || endpoints[0];
    const selCancers = () => state.cancers.filter(c => state.selectedCancers.includes(c));
    let lastSVGName = "cox_forest", lastItems = null, lastCancer = "", lastEndpoint = "OS", lastModel = "uni", lastView = "forest", lastHeatmap = null;


    container.innerHTML = `
      <div class="card">
        <div class="cx-sec">
          <h3 class="cx-h3">Genes (shared with Clinical Overview)</h3>
          <textarea id="cx-genes" class="cx-genes"></textarea>
        </div>

        <details class="cx-sec" open>
          <summary class="cx-h3" style="cursor:pointer">Cancers — 1 for a forest plot · several for an HR heatmap</summary>
          <div class="cx-toolbar" style="margin-top:8px">
            <button class="cx-mini" id="cx-c-alpha">Sort A–Z</button>
            <button class="cx-mini" id="cx-c-byn">Sort by N</button>
            <span class="cx-sep"></span>
            <button class="cx-mini" id="cx-c-all">Select all</button>
            <button class="cx-mini" id="cx-c-none">Select none</button>
          </div>
          <div id="cx-cancers" class="cx-chips"></div>
        </details>

        <div class="cx-sec">
          <div class="cx-toolbar">
            <span class="cx-ctl" style="display:inline-flex;align-items:center">Model<span class="cx-seg" id="cx-model">
              <button type="button" class="cx-seg-btn" data-model="uni">Univariate</button>
              <button type="button" class="cx-seg-btn" data-model="multi">Multivariate</button>
            </span></span>
            <label class="cx-ctl">Endpoint<select id="cx-endpoint">${endpoints.map(ep => `<option value="${ep.id}">${ep.label}</option>`).join("")}</select></label>
            <label class="cx-ctl">Gene split<select id="cx-split">
              <option value="median">Median (High/Low)</option>
              <option value="tertile">Tertile (top/bottom 1/3)</option>
              <option value="quartile">Quartile (Q4 vs Q1)</option>
            </select></label>
            <label class="cx-ctl">Months<input type="number" id="cx-months" min="0" step="1" style="width:72px" placeholder="all"></label>
          </div>
        </div>

        <div class="cx-sec">
          <div class="cx-toolbar">
            <h3 class="cx-h3" style="margin:0">Clinical factors</h3>
            <span class="cx-sep"></span>
            <button class="cx-mini" id="cx-all">Select all</button>
            <button class="cx-mini" id="cx-clear">Clear</button>
          </div>
          <div id="cx-dims" class="cx-dims"></div>
        </div>

        <div class="cx-sec">
          <div class="cx-toolbar">
            <button id="cx-run">Run</button>
            <label class="cx-ctl" id="cx-scheme-wrap" style="display:none">Colors<select id="cx-scheme"><option value="rg">Red–Green</option><option value="rb">Red–Blue</option></select></label>
            <button class="cx-mini" id="cx-swap" style="display:none">Swap rows/cols</button>
            <span class="cx-sep"></span>
            <button class="cx-mini" id="cx-csv" style="display:none">Export CSV</button>
            <button class="cx-mini" id="cx-svg" style="display:none">SVG</button>
            <button class="cx-mini" id="cx-png" style="display:none">PNG</button>
          </div>
          <div id="cx-status" class="status"></div>
          <div id="cx-result" class="cx-result"></div>
        </div>
        <div class="cx-note">Tumor only · survival endpoint selectable above. Genes are High vs Low; clinical factors are advanced vs baseline level. <span style="color:#b45309;font-weight:700">†</span> marks post-baseline factors (e.g. recurrence): off by default because they can distort baseline prognosis — enable only deliberately. Vital status is excluded (it is the survival event).</div>
      </div>`;

    const $ = s => container.querySelector(s);
    const genesEl = $("#cx-genes"), splitSel = $("#cx-split"), monthsEl = $("#cx-months"), endpointSel = $("#cx-endpoint");
    const dimsBox = $("#cx-dims"), statusEl = $("#cx-status"), resultEl = $("#cx-result"), chipBox = $("#cx-cancers");

    function commit() { saveState(state); }
    // GOIs 共享
    genesEl.value = getGOIs().join(", ");
    genesEl.addEventListener("input", () => setGOIs(parseGenes(genesEl.value)));
    const off = onGOIsChanged(list => { if (document.activeElement !== genesEl) genesEl.value = list.join(", "); });
    // 控制項初值
    endpointSel.value = state.endpoint; splitSel.value = state.split; if (state.months > 0) monthsEl.value = state.months; $("#cx-scheme").value = state.scheme;
    // Model segmented toggle（取代下拉，更明顯）
    const modelSeg = $("#cx-model");
    const setModelActive = () => modelSeg.querySelectorAll(".cx-seg-btn").forEach(b => b.classList.toggle("active", b.dataset.model === state.model));
    modelSeg.querySelectorAll(".cx-seg-btn").forEach(b => b.addEventListener("click", () => { state.model = b.dataset.model; commit(); setModelActive(); }));
    setModelActive();

    // 依目前癌種計算每個維度的可用性（base/adv 人數），無資料就灰掉並取消勾選
    function renderDims() {
      const patients = patientsInScope(dataset, selCancers());
      dimsBox.innerHTML = "";
      dims.forEach(d => {
        const r = analyzeDimension(d, patients, undefined);   // 用各級別預設指派
        const avail = r.base > 0 && r.adv > 0;
        if (!avail) state.selectedDims = state.selectedDims.filter(id => id !== d.id);
        const on = state.selectedDims.includes(d.id);
        const [bl, al] = labelsOf(d);
        const dis = DISCOURAGED[d.id];
        const chip = document.createElement("label");
        chip.className = "cx-dim" + (avail ? (on ? " on" : "") : " off") + (dis ? " warn" : "");
        chip.title = avail ? `${d.name}: ${al} vs ${bl} (n ${r.adv}/${r.base})${dis ? ` — † not recommended (${dis})` : ""}` : `${d.name}: no data in selected cancer(s)`;
        chip.innerHTML = `<input type="checkbox" ${on ? "checked" : ""} ${avail ? "" : "disabled"}><span class="cx-dim-name">${d.name}${dis ? `<sup class="cx-dag">†</sup>` : ""}</span><small>${avail ? `${al}/${bl} ${r.adv}/${r.base}` : "no data"}</small>`;
        if (avail) chip.querySelector("input").addEventListener("change", e => {
          if (e.target.checked) { if (!state.selectedDims.includes(d.id)) state.selectedDims.push(d.id); chip.classList.add("on"); }
          else { state.selectedDims = state.selectedDims.filter(id => id !== d.id); chip.classList.remove("on"); }
          commit();
        });
        dimsBox.appendChild(chip);
      });
      commit();
    }

    // 癌種 chips（多選 + 拖曳排序）
    let dragCancer = null;
    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "cx-chip" + (state.selectedCancers.includes(code) ? " on" : "");
        chip.textContent = code; chip.draggable = true;
        chip.addEventListener("click", () => {
          if (state.selectedCancers.includes(code)) state.selectedCancers = state.selectedCancers.filter(c => c !== code);
          else state.selectedCancers.push(code);
          commit(); renderCancers(); renderDims();
        });
        chip.addEventListener("dragstart", () => { dragCancer = code; chip.classList.add("dragging"); });
        chip.addEventListener("dragend", () => chip.classList.remove("dragging"));
        chip.addEventListener("dragover", e => e.preventDefault());
        chip.addEventListener("drop", e => { e.preventDefault(); state.cancers = moveInArray(state.cancers, dragCancer, code); commit(); renderCancers(); });
        chipBox.appendChild(chip);
      });
    }
    $("#cx-c-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#cx-c-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#cx-c-all").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); renderDims(); });
    $("#cx-c-none").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); renderDims(); });
    $("#cx-swap").addEventListener("click", () => { state.swap = !state.swap; commit(); if (lastView === "heatmap" && resultEl.querySelector("svg")) run(); });
    $("#cx-scheme").addEventListener("change", () => { state.scheme = $("#cx-scheme").value; commit(); if (lastView === "heatmap" && resultEl.querySelector("svg")) run(); });

    endpointSel.addEventListener("change", () => { state.endpoint = endpointSel.value; commit(); });
    splitSel.addEventListener("change", () => { state.split = splitSel.value; commit(); });
    monthsEl.addEventListener("change", () => { const v = Number(monthsEl.value); state.months = (monthsEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); });
    $("#cx-all").addEventListener("click", () => { const ps = patientsInScope(dataset, selCancers()); state.selectedDims = dims.filter(d => { if (DISCOURAGED[d.id]) return false; const r = analyzeDimension(d, ps, undefined); return r.base > 0 && r.adv > 0; }).map(d => d.id); commit(); renderDims(); });
    $("#cx-clear").addEventListener("click", () => { state.selectedDims = []; commit(); renderDims(); });

    $("#cx-run").addEventListener("click", run);
    $("#cx-csv").addEventListener("click", exportCSV);
    $("#cx-svg").addEventListener("click", downloadSVG);
    $("#cx-png").addEventListener("click", downloadPNG);

    // 取一病人的存活（套月數截斷）：回 {tm,e} 或 null
    function survOf(clin) {
      const ep = currentEndpoint();
      const os = Number(clin[ep.event]), t = Number(clin[ep.time]);
      if (!isFinite(os) || !isFinite(t) || t < 0) return null;
      const cut = state.months > 0 ? state.months : Infinity;
      let tm = t / DAYS_PER_MONTH, e = os ? 1 : 0;
      if (tm > cut) { e = 0; tm = cut; }
      return { tm, e };
    }
    // 由 (x 0/1, surv) 列表算一列 Cox 結果
    function fitRow(label, rows, tipHead) {
      const n1 = rows.filter(r => r.x === 1).length, n0 = rows.filter(r => r.x === 0).length;
      if (n1 === 0 || n0 === 0) return { label, state: "nodata", tip: `${tipHead}: one group empty (${n1}/${n0})` };
      const tm = rows.map(r => r.surv.tm), e = rows.map(r => r.surv.e), x = rows.map(r => r.x);
      const events = e.reduce((s, v) => s + v, 0);
      const cox = coxPH1(tm, e, x);
      if (!isFinite(cox.hr) || cox.hr <= 0) return { label, state: "nodata", tip: `${tipHead}: HR not estimable (n=${n1}/${n0}, events=${events})` };
      const weak = n1 < 10 || n0 < 10 || events < 10;
      return {
        label, hr: cox.hr, ciLow: cox.ciLow, ciHigh: cox.ciHigh, p: cox.p,
        n: n1 + n0, n1, n0, events, state: weak ? "weak" : "ok",
        tip: `${tipHead}: HR=${cox.hr.toFixed(2)} (${cox.ciLow.toFixed(2)}–${cox.ciHigh.toFixed(2)}), n=${n1}/${n0}, events=${events}, p=${cox.p.toPrecision(2)}`,
      };
    }

    async function run() {
      resultEl.innerHTML = ""; lastItems = null; lastHeatmap = null;
      ["cx-csv", "cx-svg", "cx-png", "cx-swap", "cx-scheme-wrap"].forEach(id => $("#" + id).style.display = "none");
      const cancers = selCancers();
      if (!cancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      const recs = [], unknown = [];
      getGOIs().forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) unknown.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      const selDims = dims.filter(d => state.selectedDims.includes(d.id));
      if (!recs.length && !selDims.length) { statusEl.className = "status err"; statusEl.textContent = "Add at least one recognized gene or select a clinical factor."; return; }

      statusEl.className = "status"; statusEl.textContent = recs.length ? `Fetching ${recs.length} gene(s)…` : "Computing…";
      let geneVals = new Map();
      try { await Promise.all(recs.map(async x => geneVals.set(x.rec.gene_id, await dataset.getGeneValues(x.rec)))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files (CORS?): " + e.message; return; }

      statusEl.textContent = "Computing…";
      const epLabel = currentEndpoint().label;
      lastEndpoint = currentEndpoint().id; lastModel = state.model;
      if (cancers.length === 1) {                                      // 單癌種 → forest
        lastView = "forest"; lastCancer = cancers[0];
        lastSVGName = `cox_${state.model}_${cancers[0]}_${lastEndpoint}`;
        const patients = patientsInScope(dataset, [cancers[0]]);
        if (state.model === "multi") renderMultivariate(recs, selDims, geneVals, patients, epLabel, unknown, cancers[0]);
        else renderUnivariate(recs, selDims, geneVals, patients, epLabel, unknown, cancers[0]);
      } else {                                                          // 多癌種 → heatmap
        lastView = "heatmap"; lastSVGName = `cox_${state.model}_heatmap_${lastEndpoint}`;
        if (state.model === "multi") drawMultiHeatmap(recs, selDims, geneVals, cancers, epLabel, unknown);
        else drawUniHeatmap(recs, selDims, geneVals, cancers, epLabel, unknown);
      }
    }

    function showExports() { $("#cx-csv").style.display = ""; $("#cx-svg").style.display = ""; $("#cx-png").style.display = ""; }

    // ── Univariate：每個因子各自一個 Cox ──
    function renderUnivariate(recs, selDims, geneVals, patients, epLabel, unknown, cancer) {
      const items = [];
      recs.forEach(g => {
        const vals = geneVals.get(g.rec.gene_id);
        const elig = patients.map(p => ({ p, expr: vals[p.idx], surv: survOf(p.clin) })).filter(o => isFinite(o.expr) && o.surv);
        const grp = splitGroups(elig.map(o => o.expr), state.split);
        const rows = [];
        elig.forEach((o, i) => { if (grp[i] === 1) rows.push({ x: 1, surv: o.surv }); else if (grp[i] === 0) rows.push({ x: 0, surv: o.surv }); });
        items.push(fitRow(`${g.label}: High vs Low`, rows, g.label));
      });
      selDims.forEach(d => {
        const [bl, al] = labelsOf(d);
        const rows = [];
        patients.forEach(p => {
          const surv = survOf(p.clin); if (!surv) return;
          const b = classify(d, p.clin[d.field], undefined);
          if (b === "advanced") rows.push({ x: 1, surv }); else if (b === "baseline") rows.push({ x: 0, surv });
        });
        items.push(fitRow(`${d.name}: ${al} vs ${bl}`, rows, d.name));
      });
      // FDR：跨可估計的列校正
      const est = items.filter(it => it.state !== "nodata" && isFinite(it.p));
      const q = benjaminiHochberg(est.map(it => it.p));
      est.forEach((it, k) => { it.q = q[k]; const st = pStars(q[k]); it.stars = st === "ns" ? "" : st; it.pText = fmtP(it.p); it.tip += `, q=${q[k].toPrecision(2)}`; });
      items.forEach(it => { if (it.state === "nodata") it.pText = ""; });

      lastItems = items;
      resultEl.innerHTML = forestSVG(items, {
        caption: `Univariate Cox · ${cancer} · ${epLabel}${state.months > 0 ? ` · ${state.months}-mo` : ""}`,
        scheme: "rg",
      });
      const note = document.createElement("div"); note.className = "cx-legend";
      note.textContent = `Each factor in its own univariate Cox (${epLabel}). Genes: High vs Low (${state.split}); clinical: advanced vs baseline. ★ FDR q<0.05 ★★<0.01 ★★★<0.001 ★★★★<0.0001 · faded = <10 per group or <10 events · red HR>1 (worse), green HR<1 (better)${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      showExports();
      statusEl.textContent = `Done — ${items.length} factor(s).`;
    }

    // ── Multivariate：所有選到的因子放進同一個 Cox（互相校正，complete-case）──
    function renderMultivariate(recs, selDims, geneVals, patients, epLabel, unknown, cancer) {
      const elig = patients.map(p => ({ p, surv: survOf(p.clin) })).filter(o => o.surv);
      const covs = [];   // 每個共變量：{label, code:[每位 elig 病人的 0/1 或 null=缺]}
      recs.forEach(g => {
        const vals = geneVals.get(g.rec.gene_id);
        const expr = elig.map(o => vals[o.p.idx]);
        const finiteIdx = []; expr.forEach((v, i) => { if (isFinite(v)) finiteIdx.push(i); });
        const grpF = splitGroups(finiteIdx.map(i => expr[i]), state.split);   // split 在「有表現」的病人上算
        const code = new Array(elig.length).fill(null);
        finiteIdx.forEach((i, k) => { code[i] = grpF[k] === 1 ? 1 : grpF[k] === 0 ? 0 : null; });
        covs.push({ label: `${g.label}: High vs Low`, code });
      });
      selDims.forEach(d => {
        const [bl, al] = labelsOf(d);
        const code = elig.map(o => { const b = classify(d, o.p.clin[d.field], undefined); return b === "advanced" ? 1 : b === "baseline" ? 0 : null; });
        covs.push({ label: `${d.name}: ${al} vs ${bl}`, code });
      });
      if (!covs.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one factor for the model."; return; }
      const keep = [];
      for (let i = 0; i < elig.length; i++) if (covs.every(c => c.code[i] != null)) keep.push(i);
      const nEvents = keep.reduce((s, i) => s + elig[i].surv.e, 0);
      if (keep.length < covs.length + 2 || nEvents < 2) {
        statusEl.className = "status err";
        statusEl.textContent = `Not enough complete cases (n=${keep.length}, events=${nEvents}) for ${covs.length} covariate(s). Remove a factor.`;
        return;
      }
      const times = keep.map(i => elig[i].surv.tm), events = keep.map(i => elig[i].surv.e);
      const X = keep.map(i => covs.map(c => c.code[i]));
      const fit = coxPH(times, events, X);
      if (fit.error) { statusEl.className = "status err"; statusEl.textContent = `Model not estimable: ${fit.error}. Try removing a collinear factor.`; return; }

      const items = covs.map((c, k) => {
        const hr = fit.hr[k], lo = fit.ciLow[k], hi = fit.ciHigh[k], p = fit.p[k];
        const st = pStars(p);
        return {
          label: c.label, hr, ciLow: lo, ciHigh: hi, p, pText: fmtP(p), stars: st === "ns" ? "" : st, state: "ok",
          tip: `${c.label}: adjusted HR=${hr.toFixed(2)} (${lo.toFixed(2)}–${hi.toFixed(2)}), p=${p.toPrecision(2)}`,
        };
      });
      lastItems = items;
      resultEl.innerHTML = forestSVG(items, {
        caption: `Multivariate Cox · ${cancer} · ${epLabel}${state.months > 0 ? ` · ${state.months}-mo` : ""}`,
        scheme: "rg",
      });
      const epv = nEvents / covs.length;
      const note = document.createElement("div"); note.className = "cx-legend";
      note.innerHTML = `One multivariate Cox (${epLabel}), all factors mutually adjusted · n=${fit.n} complete cases, ${fit.events} events, ${covs.length} covariates. Genes High vs Low (${state.split}); clinical advanced vs baseline. p = Wald. ★ p<0.05 ★★<0.01 ★★★<0.001 ★★★★<0.0001 · red adjusted HR>1 (worse), green <1 (better)${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`
        + (epv < 10 ? `<br><span style="color:#b45309">${epv.toFixed(1)} events per covariate (&lt;10) — estimates may be unstable; consider fewer factors.</span>` : "");
      resultEl.appendChild(note);
      showExports();
      statusEl.textContent = `Done — multivariate, n=${fit.n}, ${fit.events} events.`;
    }

    // ── Univariate × 多癌種 → HR heatmap（列=factor、欄=癌種、色=log2 HR；FDR 每癌種欄內）──
    function drawUniHeatmap(recs, selDims, geneVals, cancers, epLabel, unknown) {
      const factors = [
        ...recs.map(g => ({ kind: "gene", g, label: `${g.label}: High vs Low` })),
        ...selDims.map(d => { const [bl, al] = labelsOf(d); return { kind: "dim", d, label: `${d.name}: ${al} vs ${bl}` }; }),
      ];
      const byCancer = {}; cancers.forEach(c => byCancer[c] = patientsInScope(dataset, [c]));
      const cells = [];                 // cells[factorIdx][cancerIdx]
      const pPerCol = cancers.map(() => []);
      factors.forEach((f) => {
        const rowCells = [];
        cancers.forEach((c, ci) => {
          const patients = byCancer[c], rows = [];
          if (f.kind === "gene") {
            const vals = geneVals.get(f.g.rec.gene_id);
            const elig = patients.map(p => ({ p, expr: vals[p.idx], surv: survOf(p.clin) })).filter(o => isFinite(o.expr) && o.surv);
            const grp = splitGroups(elig.map(o => o.expr), state.split);
            elig.forEach((o, i) => { if (grp[i] === 1) rows.push({ x: 1, surv: o.surv }); else if (grp[i] === 0) rows.push({ x: 0, surv: o.surv }); });
          } else {
            patients.forEach(p => { const surv = survOf(p.clin); if (!surv) return; const b = classify(f.d, p.clin[f.d.field], undefined); if (b === "advanced") rows.push({ x: 1, surv }); else if (b === "baseline") rows.push({ x: 0, surv }); });
          }
          const fit = fitRow(f.label, rows, `${f.label} / ${c}`);
          if (fit.state === "nodata") rowCells.push({ state: "nodata", tip: fit.tip });
          else { rowCells.push({ value: Math.log2(fit.hr), state: fit.state, stars: "", tip: fit.tip, _hr: fit.hr }); pPerCol[ci].push({ ri: cells.length, p: fit.p }); }
        });
        cells.push(rowCells);
      });
      cancers.forEach((c, ci) => {   // FDR：每癌種欄內各自校正
        const q = benjaminiHochberg(pPerCol[ci].map(x => x.p));
        pPerCol[ci].forEach((x, k) => { const cell = cells[x.ri][ci]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      });
      let maxAbs = 0; cells.forEach(row => row.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) maxAbs = Math.max(maxAbs, Math.abs(c.value)); }));
      const colorMax = Math.min(4, Math.max(0.5, Math.ceil(maxAbs * 10) / 10));
      lastHeatmap = { factors, cancers, cells };
      const facLab = factors.map(f => f.label);
      let rLab = facLab, cLab = cancers, grid = cells;
      if (state.swap) { rLab = cancers; cLab = facLab; grid = cancers.map((_, ci) => factors.map((_, ri) => cells[ri][ci])); }
      resultEl.innerHTML = heatmapSVG(rLab, cLab, grid, { colorMax, scheme: state.scheme, legendLabel: "log2 HR", caption: `Univariate Cox HR · ${epLabel}${state.months > 0 ? ` · ${state.months}-mo` : ""}` });
      const note = document.createElement("div"); note.className = "cx-legend";
      note.textContent = `Per-factor univariate Cox HR per cancer (${epLabel}). red HR>1 (worse), ${state.scheme === "rb" ? "blue" : "green"} HR<1 (better) · ★ FDR within each cancer · faded = <10 per group or <10 events · grey = not estimable${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      ["cx-csv", "cx-svg", "cx-png", "cx-swap", "cx-scheme-wrap"].forEach(id => $("#" + id).style.display = "");
      statusEl.textContent = `Done — ${factors.length} factor(s) × ${cancers.length} cancer(s).`;
    }

    // 某癌種的多變量模型 → 每 factor 的 adjusted HR（該癌種無資料/單組的 factor 不納入模型）
    function cancerMultiHR(cancer, factors, geneVals) {
      const patients = patientsInScope(dataset, [cancer]);
      const elig = patients.map(p => ({ p, surv: survOf(p.clin) })).filter(o => o.surv);
      const used = [];
      factors.forEach((f, fi) => {
        let code;
        if (f.kind === "gene") {
          const vals = geneVals.get(f.g.rec.gene_id);
          const expr = elig.map(o => vals[o.p.idx]);
          const finiteIdx = []; expr.forEach((v, i) => { if (isFinite(v)) finiteIdx.push(i); });
          const grpF = splitGroups(finiteIdx.map(i => expr[i]), state.split);
          code = new Array(elig.length).fill(null);
          finiteIdx.forEach((i, k) => { code[i] = grpF[k] === 1 ? 1 : grpF[k] === 0 ? 0 : null; });
        } else {
          code = elig.map(o => { const b = classify(f.d, o.p.clin[f.d.field], undefined); return b === "advanced" ? 1 : b === "baseline" ? 0 : null; });
        }
        const n1 = code.reduce((s, c) => s + (c === 1 ? 1 : 0), 0), n0 = code.reduce((s, c) => s + (c === 0 ? 1 : 0), 0);
        if (n1 > 0 && n0 > 0) used.push({ fi, code });
      });
      if (used.length < 1) return { byFactor: {} };
      const keep = [];
      for (let i = 0; i < elig.length; i++) if (used.every(u => u.code[i] != null)) keep.push(i);
      const nEvents = keep.reduce((s, i) => s + elig[i].surv.e, 0);
      if (keep.length < used.length + 2 || nEvents < 2) return { byFactor: {} };
      const times = keep.map(i => elig[i].surv.tm), events = keep.map(i => elig[i].surv.e);
      const X = keep.map(i => used.map(u => u.code[i]));
      const fit = coxPH(times, events, X);
      if (fit.error) return { byFactor: {} };
      const epv = nEvents / used.length, byFactor = {};
      used.forEach((u, k) => { byFactor[u.fi] = { hr: fit.hr[k], p: fit.p[k], weak: epv < 10 }; });
      return { byFactor };
    }

    // ── Multivariate × 多癌種 → adjusted HR heatmap（每癌種各建一個模型）──
    function drawMultiHeatmap(recs, selDims, geneVals, cancers, epLabel, unknown) {
      const factors = [
        ...recs.map(g => ({ kind: "gene", g, label: `${g.label}: High vs Low` })),
        ...selDims.map(d => { const [bl, al] = labelsOf(d); return { kind: "dim", d, label: `${d.name}: ${al} vs ${bl}` }; }),
      ];
      const perCancer = cancers.map(c => cancerMultiHR(c, factors, geneVals));
      const cells = factors.map((f, fi) => cancers.map((c, ci) => {
        const r = perCancer[ci].byFactor[fi];
        if (!r || !isFinite(r.hr) || r.hr <= 0) return { state: "nodata", tip: `${f.label} / ${c}: not in model / not estimable` };
        return { value: Math.log2(r.hr), state: r.weak ? "weak" : "ok", stars: "", tip: `${f.label} / ${c}: adj HR=${r.hr.toFixed(2)}, p=${r.p.toPrecision(2)}`, _hr: r.hr, _p: r.p };
      }));
      cancers.forEach((c, ci) => {   // FDR：每癌種欄內各自校正（Wald p）
        const idxs = []; factors.forEach((f, fi) => { const cell = cells[fi][ci]; if (cell.state !== "nodata") idxs.push({ fi, p: cell._p }); });
        const q = benjaminiHochberg(idxs.map(x => x.p));
        idxs.forEach((x, k) => { const cell = cells[x.fi][ci]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      });
      let maxAbs = 0; cells.forEach(row => row.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) maxAbs = Math.max(maxAbs, Math.abs(c.value)); }));
      const colorMax = Math.min(4, Math.max(0.5, Math.ceil(maxAbs * 10) / 10));
      lastHeatmap = { factors, cancers, cells };
      const facLab = factors.map(f => f.label);
      let rLab = facLab, cLab = cancers, grid = cells;
      if (state.swap) { rLab = cancers; cLab = facLab; grid = cancers.map((_, ci) => factors.map((_, ri) => cells[ri][ci])); }
      resultEl.innerHTML = heatmapSVG(rLab, cLab, grid, { colorMax, scheme: state.scheme, legendLabel: "log2 adj HR", caption: `Multivariate Cox adjusted HR · ${epLabel}${state.months > 0 ? ` · ${state.months}-mo` : ""}` });
      const note = document.createElement("div"); note.className = "cx-legend";
      note.textContent = `Each cancer = its own multivariate Cox (${epLabel}), adjusted HR per factor. Exploratory pan-cancer — each model adjusts for the factors available in that cancer (sets may differ), so columns are not strictly comparable. red adj.HR>1, ${state.scheme === "rb" ? "blue" : "green"} <1 · ★ FDR within each cancer · faded = EPV<10 · grey = not in model${unknown.length ? ` · Unrecognized: ${unknown.join(", ")}` : ""}`;
      resultEl.appendChild(note);
      ["cx-csv", "cx-svg", "cx-png", "cx-swap", "cx-scheme-wrap"].forEach(id => $("#" + id).style.display = "");
      statusEl.textContent = `Done — multivariate heatmap, ${factors.length} factor(s) × ${cancers.length} cancer(s).`;
    }

    function exportCSV() {
      if (lastView === "heatmap") {
        if (!lastHeatmap) return;
        const { factors, cancers, cells } = lastHeatmap;
        let csv = "Factor," + cancers.map(c => `"${c} HR"`).join(",") + "\n";
        factors.forEach((f, ri) => { csv += `"${f.label}",` + cancers.map((c, ci) => { const cell = cells[ri][ci]; return cell && cell.state !== "nodata" && cell._hr != null ? cell._hr.toFixed(4) : "NA"; }).join(",") + "\n"; });
        dl(new Blob([csv], { type: "text/csv" }), `${lastSVGName}.csv`);
        return;
      }
      if (!lastItems) return;
      let csv, name;
      if (lastModel === "multi") {
        csv = "Factor,adjHR,CI_low,CI_high,p_wald\n";
        lastItems.forEach(it => { csv += `"${it.label}",${it.hr.toFixed(4)},${it.ciLow.toFixed(4)},${it.ciHigh.toFixed(4)},${it.p.toExponential(4)}\n`; });
        name = `cox_multivariate_${lastCancer}_${lastEndpoint}.csv`;
      } else {
        csv = "Factor,HR,CI_low,CI_high,p,q,n_index,n_ref,events\n";
        lastItems.forEach(it => {
          if (it.state === "nodata") { csv += `"${it.label}",NA,NA,NA,NA,NA,NA,NA,NA\n`; return; }
          csv += `"${it.label}",${it.hr.toFixed(4)},${it.ciLow.toFixed(4)},${it.ciHigh.toFixed(4)},${it.p.toExponential(4)},${(it.q != null ? it.q.toExponential(4) : "NA")},${it.n1},${it.n0},${it.events}\n`;
        });
        name = `cox_univariate_${lastCancer}_${lastEndpoint}.csv`;
      }
      dl(new Blob([csv], { type: "text/csv" }), name);
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
    renderDims();
  },
};
