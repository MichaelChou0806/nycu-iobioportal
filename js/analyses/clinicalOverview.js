// =====================================================================
// analyses/clinicalOverview.js —「Clinical Overview」（交付 A + B）
// ---------------------------------------------------------------------
// A：多 GOI / 癌種多選 / 維度勾選 + 逐級別指派 / 即時人數 + 黃旗 /
//    拖曳排序 + 快捷排序 / 自動記憶 / 匯出入設定 / 命名清單
// B：Run -> 對每個 (基因×癌種×維度) 跑 U-test，畫成：
//    - Expanded：列=基因、欄=癌種，色=log2FC，★=該癌種內檢定(FDR)
//    - Condensed：所有癌種合併，列=基因、欄=維度兩側，色=row z-score，★在 advanced 側
//    缺資料三態（grey 斜線 / 淡化小點 / 正常）；最終圖無警告 flag；可匯出結果表
// =====================================================================

import { loadDimensions, patientsInScope, analyzeDimension, warnFor, classify }
  from "../core/dimensions.js";
import { loadLast, saveLast, exportState, importStateFromFile,
         saveNamed, loadNamed, listNames, deleteNamed, reconcile }
  from "../core/state.js";
import { mannWhitney, median, mean, sd, log2FC, benjaminiHochberg, zscoreRow, pStars }
  from "../core/stats.js";
import { heatmapSVG, multiBarSVG } from "../core/plots.js";
import { getGOIs, setGOIs, onGOIsChanged, parseGenes } from "../core/gois.js";

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .co-sec{margin-bottom:18px}
    .co-sec h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .co-sum{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;font-weight:600;cursor:pointer;padding:4px 0}
    details.co-sec>summary{margin-bottom:6px}
    .co-genes{width:100%;min-height:56px;padding:9px 11px;border:1px solid var(--line);border-radius:7px;font:13px/1.5 monospace;resize:vertical}
    .co-toolbar{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:8px;align-items:center}
    .co-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .co-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .co-chips{display:flex;flex-wrap:wrap;gap:6px}
    .co-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .co-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .co-chip.dragging,.co-dim.dragging{opacity:.4}
    .co-dim{display:flex;flex-wrap:wrap;align-items:center;gap:10px;padding:8px 0;border-bottom:1px solid #f0f2f5}
    .co-dim.off{opacity:.45}
    .co-handle{cursor:grab;color:#9aa3ad;font-size:16px;user-select:none;padding:0 2px}
    .co-dim .name{font-weight:600;min-width:84px}
    .co-cov{font-size:12px;color:var(--muted)}
    .co-cnt{font-size:12px;font-weight:600}
    .co-warn{font-size:11.5px;background:#fffbeb;color:#92400e;border:1px solid #fde68a;border-radius:6px;padding:2px 8px}
    .co-maperr{font-size:11.5px;background:#fef2f2;color:#991b1b;border:1px solid #fecaca;border-radius:6px;padding:2px 8px}
    .co-levels{display:flex;flex-wrap:wrap;gap:8px;width:100%;padding:4px 0 0 110px}
    .co-levels .lv{display:flex;align-items:center;gap:4px;font-size:12px}
    .co-levels select{min-width:auto;padding:3px 6px;font-size:12px}
    .co-io{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:4px;padding:10px 12px;background:#f8fafc;border:1px solid var(--line);border-radius:8px}
    .co-io input[type=text],.co-io select{min-width:120px;padding:5px 8px;font-size:12px}
    .co-note{font-size:11.5px;color:var(--muted)}
    .co-result{overflow-x:auto}
    .co-result svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .co-disabled{color:#9aa3ad;font-style:italic}
    .co-legend{font-size:11.5px;color:var(--muted);margin-top:6px}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

function moveInArray(arr, fromId, toId) {
  if (fromId == null || fromId === toId) return arr.slice();
  const a = arr.slice(); const fi = a.indexOf(fromId); if (fi < 0) return a;
  a.splice(fi, 1); let ti = a.indexOf(toId); if (ti < 0) ti = a.length; a.splice(ti, 0, fromId); return a;
}
// 取維度兩側的顯示標籤（自明）：
//  - numericSplit（Age）：依目前 cutoff → "<60" / ">=60"
//  - ordinal（T/Stage/Grade）：依目前每級別歸屬 → 範圍標籤，如 "T1–T3" / "T4"；不連續用逗號 "T1,T2,T4"
//  - 其他 binary：用定義檔的 label（已改為自明，如 Alcohol- / Alcohol+）
function labelsOf(dim, assignment) {
  if (dim.type === "binary") {
    if (dim.numericSplit) {
      const c = (assignment && assignment.cutoff != null) ? assignment.cutoff : dim.numericSplit.cutoff;
      return [`<${c}`, `>=${c}`];
    }
    if (dim.baseline && dim.advanced) return [dim.baseline.label, dim.advanced.label];
  }
  if (dim.type === "ordinal") {
    return [ordinalRange(dim, assignment, "baseline"), ordinalRange(dim, assignment, "advanced")];
  }
  return ["baseline", "advanced"];
}
function ordinalRange(dim, assignment, side) {
  const prefix = dim.labelPrefix || "";
  const order = dim.levels.map(l => l.id);
  const ids = dim.levels.filter(lv => ((assignment && assignment[lv.id]) || lv.default) === side).map(lv => lv.id);
  if (!ids.length) return side === "baseline" ? "Baseline" : "Advanced";
  const idxs = ids.map(id => order.indexOf(id)).sort((a, b) => a - b);
  const contiguous = idxs.every((v, i) => i === 0 || v === idxs[i - 1] + 1);
  const range = (contiguous && ids.length > 1) ? `${ids[0]}–${ids[ids.length - 1]}` : ids.join(",");
  return prefix + range;
}
// 只顯示顯著星號；不顯著回空字串（不標 ns）
function starsOf(q) { const s = pStars(q); return (s === "ns") ? "" : s; }

export const clinicalOverview = {
  id: "clinicalOverview",
  name: "Clinical Overview",

  async mount(container, { dataset }) {
    injectStyles();
    container.innerHTML = `<div class="card"><div class="status">Loading dimension definitions…</div></div>`;

    const dimUrl = (dataset.config && dataset.config.dimensionsUrl) || "config/dimensions.tcga.json";
    let DEF;
    try { DEF = await loadDimensions(dimUrl); }
    catch (e) { container.innerHTML = `<div class="card"><div class="status err">Failed to load dimension definitions: ${e.message}</div></div>`; return; }
    const THR = DEF.warnThresholds || { minPerGroup: 10, maxRatio: 10 };
    const dimById = Object.fromEntries(DEF.dimensions.map(d => [d.id, d]));
    const avail = { cancers: dataset.cancers.map(c => c.code), dimIds: DEF.dimensions.map(d => d.id) };
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));

    let state = reconcile(loadLast(), avail);
    const _shared = getGOIs();
    if (_shared.length) state.genes = _shared; else if (state.genes.length) setGOIs(state.genes);
    let dragCancer = null, dragDim = null;
    let lastResults = null;

    container.innerHTML = `
      <div class="card">
        <div class="co-sec">
          <h3>Genes (comma / space / newline separated)</h3>
          <textarea id="co-genes" class="co-genes"></textarea>
        </div>

        <details class="co-sec" open>
          <summary class="co-sum">Cancers — drag to reorder</summary>
          <div class="co-toolbar" style="margin-top:8px">
            <button class="co-mini" id="cs-alpha">Sort A–Z</button>
            <button class="co-mini" id="cs-byn">Sort by N</button>
            <button class="co-mini" id="cs-def">Default order</button>
            <span class="co-sep"></span>
            <button class="co-mini" id="cs-all">Select all</button>
            <button class="co-mini" id="cs-none">Select none</button>
          </div>
          <div id="co-cancers" class="co-chips"></div>
        </details>

        <details class="co-sec" open>
          <summary class="co-sum">Clinical dimensions — drag (≡) to reorder; tick to include</summary>
          <div class="co-toolbar" style="margin-top:8px"><button class="co-mini" id="ds-all">Select all</button><button class="co-mini" id="ds-none">Select none</button><span class="co-sep"></span><button class="co-mini" id="ds-def">Default order</button></div>
          <div id="co-dims"></div>
        </details>

        <div class="co-sec">
          <div class="co-toolbar">
            <label class="co-note">View
              <select id="co-mode">
                <option value="expanded">Expanded · genes × cancers (one dimension)</option>
                <option value="condensed">Condensed · genes × dimensions (cancers merged)</option>
              </select>
            </label>
            <label class="co-note" id="co-dimpick-wrap">Dimension
              <select id="co-dimpick"></select>
            </label>
            <button id="co-run">Run</button>
            <button class="co-mini" id="co-expres" style="display:none">Export results (CSV)</button>
          </div>
          <div id="co-status" class="status"></div>
          <div id="co-result" class="co-result"></div>
        </div>

        <details style="margin-top:6px">
          <summary class="co-note" style="cursor:pointer">Settings &amp; saved lists</summary>
          <div class="co-io" style="margin-top:8px">
            <button class="co-mini" id="co-export">Export</button>
            <label class="co-mini" style="display:inline-block">Import<input type="file" id="co-import" accept=".json" style="display:none"></label>
            <span class="co-sep"></span>
            <input type="text" id="co-savename" placeholder="Save current as…">
            <button class="co-mini" id="co-save">Save</button>
            <select id="co-loadsel"><option value="">— saved —</option></select>
            <button class="co-mini" id="co-load">Load</button>
            <button class="co-mini" id="co-del">Delete</button>
          </div>
        </details>
      </div>`;

    const $ = s => container.querySelector(s);
    const chipBox = $("#co-cancers"), dimsBox = $("#co-dims"), genesEl = $("#co-genes");
    const modeSel = $("#co-mode"), dimPick = $("#co-dimpick"), dimPickWrap = $("#co-dimpick-wrap");
    const statusEl = $("#co-status"), resultEl = $("#co-result");

    function commit() { saveLast(state); }
    function syncGenes() { genesEl.value = state.genes.join(", "); }
    function selectedCancersOrdered() { return state.cancers.filter(c => state.selectedCancers.includes(c)); }
    function selectedDimsOrdered() { return state.dimensions.filter(id => state.selectedDims.includes(id)).map(id => dimById[id]).filter(Boolean); }
    function assignmentOf(dim) {
      if (dim.type === "ordinal") {
        const a = {}; dim.levels.forEach(lv => { a[lv.id] = (state.ordinalAssign[dim.id] && state.ordinalAssign[dim.id][lv.id]) || lv.default; }); return a;
      }
      if (dim.numericSplit) {
        return { cutoff: (state.numericCutoff && state.numericCutoff[dim.id] != null) ? state.numericCutoff[dim.id] : dim.numericSplit.cutoff };
      }
      return null;
    }

    // ---------- 渲染：癌種 ----------
    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "co-chip" + (state.selectedCancers.includes(code) ? " on" : "");
        chip.dataset.code = code; chip.textContent = code; chip.draggable = true;
        chip.addEventListener("click", () => {
          if (state.selectedCancers.includes(code)) { state.selectedCancers = state.selectedCancers.filter(c => c !== code); chip.classList.remove("on"); }
          else { state.selectedCancers.push(code); chip.classList.add("on"); }
          commit(); recompute();
        });
        chip.addEventListener("dragstart", () => { dragCancer = code; chip.classList.add("dragging"); });
        chip.addEventListener("dragend", () => chip.classList.remove("dragging"));
        chip.addEventListener("dragover", e => e.preventDefault());
        chip.addEventListener("drop", e => { e.preventDefault(); state.cancers = moveInArray(state.cancers, dragCancer, code); commit(); renderCancers(); });
        chipBox.appendChild(chip);
      });
    }

    // ---------- 渲染：維度 ----------
    function renderDims() {
      dimsBox.innerHTML = "";
      state.dimensions.forEach(id => {
        const dim = dimById[id]; if (!dim) return;
        const row = document.createElement("div"); row.className = "co-dim"; row.dataset.id = id;
        let levelsHtml = "";
        if (dim.type === "ordinal") {
          const a = assignmentOf(dim);
          levelsHtml = `<div class="co-levels">` + dim.levels.map(lv => `
            <span class="lv">${lv.label}
              <select data-dim="${id}" data-level="${lv.id}">
                <option value="baseline"${a[lv.id] === "baseline" ? " selected" : ""}>Baseline</option>
                <option value="advanced"${a[lv.id] === "advanced" ? " selected" : ""}>Advanced</option>
                <option value="ignore"${a[lv.id] === "ignore" ? " selected" : ""}>Ignore</option>
              </select></span>`).join("") + `</div>`;
        } else if (dim.numericSplit) {
          const c = assignmentOf(dim).cutoff;
          levelsHtml = `<div class="co-levels"><span class="lv">cutoff <input type="number" class="co-cut" value="${c}" style="width:64px;padding:3px 6px;font-size:12px"></span></div>`;
        }
        row.innerHTML = `
          <span class="co-handle" draggable="true" title="拖曳排序">≡</span>
          <label style="display:flex;align-items:center;gap:6px">
            <input type="checkbox" class="co-on" ${state.selectedDims.includes(id) ? "checked" : ""}> <span class="name">${dim.name}</span>
          </label>
          <span class="co-cov" id="cov-${id}"></span>
          <span class="co-cnt" id="cnt-${id}"></span>
          <span class="co-warn hidden" id="warn-${id}"></span>
          <span class="co-maperr hidden" id="map-${id}"></span>
          ${levelsHtml}`;
        dimsBox.appendChild(row);

        row.querySelector(".co-on").addEventListener("change", e => {
          if (e.target.checked) { if (!state.selectedDims.includes(id)) state.selectedDims.push(id); }
          else state.selectedDims = state.selectedDims.filter(x => x !== id);
          commit(); refreshDimPick();
        });
        if (dim.type === "ordinal") {
          row.querySelectorAll("select").forEach(sel => sel.addEventListener("change", () => {
            if (!state.ordinalAssign[id]) state.ordinalAssign[id] = {};
            state.ordinalAssign[id][sel.dataset.level] = sel.value; commit(); recompute();
          }));
        }
        if (dim.numericSplit) {
          const inp = row.querySelector(".co-cut");
          if (inp) inp.addEventListener("change", () => {
            const v = Number(inp.value);
            if (isFinite(v)) { if (!state.numericCutoff) state.numericCutoff = {}; state.numericCutoff[dim.id] = v; commit(); recompute(); }
          });
        }
        const handle = row.querySelector(".co-handle");
        handle.addEventListener("dragstart", () => { dragDim = id; row.classList.add("dragging"); });
        handle.addEventListener("dragend", () => row.classList.remove("dragging"));
        row.addEventListener("dragover", e => e.preventDefault());
        row.addEventListener("drop", e => { e.preventDefault(); state.dimensions = moveInArray(state.dimensions, dragDim, id); commit(); renderDims(); recompute(); refreshDimPick(); });
      });
    }

    // ---------- 即時覆蓋/人數/黃旗 ----------
    function recompute() {
      const patients = state.selectedCancers.length ? patientsInScope(dataset, state.selectedCancers) : [];
      DEF.dimensions.forEach(dim => {
        const row = dimsBox.querySelector(`.co-dim[data-id="${dim.id}"]`); if (!row) return;
        const covEl = $(`#cov-${dim.id}`), cntEl = $(`#cnt-${dim.id}`), warnEl = $(`#warn-${dim.id}`), mapEl = $(`#map-${dim.id}`);
        const onBox = row.querySelector(".co-on");
        if (!patients.length) { covEl.textContent = "—"; cntEl.textContent = ""; warnEl.classList.add("hidden"); mapEl.classList.add("hidden"); return; }
        const r = analyzeDimension(dim, patients, assignmentOf(dim));
        if (r.nField === 0) {
          row.classList.add("off"); onBox.checked = false; onBox.disabled = true;
          state.selectedDims = state.selectedDims.filter(x => x !== dim.id);
          covEl.innerHTML = `<span class="co-disabled">no data in selection</span>`;
          cntEl.textContent = ""; warnEl.classList.add("hidden"); mapEl.classList.add("hidden"); return;
        }
        row.classList.remove("off"); onBox.disabled = false;
        covEl.textContent = `${r.nCancers} cancer(s) · ${r.nField} pts`;
        const lbl = labelsOf(dim, assignmentOf(dim));
        cntEl.textContent = `${lbl[0]}=${r.base} / ${lbl[1]}=${r.adv}`;
        if (r.nMapped < r.nField * 0.8) { mapEl.classList.remove("hidden"); mapEl.textContent = `only ${r.nMapped}/${r.nField} mapped — check definition`; }
        else mapEl.classList.add("hidden");
        const w = warnFor(r.base, r.adv, THR);
        if (w.length) { warnEl.classList.remove("hidden"); warnEl.textContent = "⚠ " + w.join("; "); } else warnEl.classList.add("hidden");
      });
      refreshDimPick();
    }

    // 認得的基因數（決定單基因/多基因的呈現方式）
    function recognizedGeneCount() {
      let n = 0; state.genes.forEach(g => { const r = dataset.resolveGene(g); if (!r.error && !r.multiple) n++; }); return n;
    }
    // 展開模式的維度下拉：單基因時隱藏（單基因改看所有維度），多基因才需選一個維度
    function refreshDimPick() {
      const dims = selectedDimsOrdered();
      const cur = dimPick.value;
      dimPick.innerHTML = dims.map(d => `<option value="${d.id}">${d.name}</option>`).join("");
      if (dims.some(d => d.id === cur)) dimPick.value = cur;
      const single = recognizedGeneCount() === 1;
      dimPickWrap.style.display = (state.mode === "expanded" && !single) ? "" : "none";
      const expOpt = modeSel.querySelector('option[value="expanded"]');
      if (expOpt) expOpt.textContent = single ? "Expanded · dimensions × cancers" : "Expanded · genes × cancers (one dimension)";
    }

    function renderAll() { renderCancers(); renderDims(); recompute(); }

    // ---------- 快捷排序/選取 ----------
    $("#cs-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#cs-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#cs-def").addEventListener("click", () => { state.cancers = avail.cancers.slice(); commit(); renderCancers(); });
    $("#cs-all").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); recompute(); });
    $("#cs-none").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); recompute(); });
    $("#ds-all").addEventListener("click", () => { state.selectedDims = state.dimensions.slice(); commit(); renderDims(); recompute(); refreshDimPick(); });
    $("#ds-none").addEventListener("click", () => { state.selectedDims = []; commit(); renderDims(); recompute(); refreshDimPick(); });
    $("#ds-def").addEventListener("click", () => { state.dimensions = avail.dimIds.slice(); commit(); renderDims(); recompute(); });

    genesEl.addEventListener("input", () => { state.genes = parseGenes(genesEl.value); setGOIs(state.genes); commit(); refreshDimPick(); });
    modeSel.value = state.mode || "expanded";
    modeSel.addEventListener("change", () => { state.mode = modeSel.value; commit(); refreshDimPick(); });

    // ---------- 設定匯出入 / 命名清單 ----------
    function refreshNamed() { const sel = $("#co-loadsel"); const cur = sel.value; sel.innerHTML = `<option value="">— saved —</option>` + listNames().map(n => `<option value="${n}">${n}</option>`).join(""); sel.value = cur; }
    function applyState(s) { state = reconcile(s, avail); syncGenes(); modeSel.value = state.mode || "expanded"; renderAll(); }
    $("#co-export").addEventListener("click", () => exportState(state, "tcga-overview-settings.json"));
    $("#co-import").addEventListener("change", async e => { if (!e.target.files[0]) return; try { applyState(await importStateFromFile(e.target.files[0])); commit(); } catch (err) { alert(err.message); } e.target.value = ""; });
    $("#co-save").addEventListener("click", () => { const n = $("#co-savename").value.trim(); if (!n) { alert("Enter a list name first."); return; } saveNamed(n, state); refreshNamed(); $("#co-loadsel").value = n; });
    $("#co-load").addEventListener("click", () => { const n = $("#co-loadsel").value; if (!n) return; const s = loadNamed(n); if (s) { applyState(s); commit(); } });
    $("#co-del").addEventListener("click", () => { const n = $("#co-loadsel").value; if (!n) return; if (confirm(`Delete list "${n}"?`)) { deleteNamed(n); refreshNamed(); } });

    // =====================================================================
    // RUN：計算 + 繪圖
    // =====================================================================
    $("#co-run").addEventListener("click", run);
    $("#co-expres").addEventListener("click", exportResults);

    async function run() {
      resultEl.innerHTML = ""; $("#co-expres").style.display = "none"; lastResults = null;
      const cancers = selectedCancersOrdered();
      const dims = selectedDimsOrdered();
      if (!cancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      if (!dims.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one clinical dimension."; return; }

      // 基因解析
      const recs = []; const unknown = [];
      state.genes.forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) unknown.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (!recs.length) { statusEl.className = "status err"; statusEl.textContent = "No recognized genes." + (unknown.length ? ` Unrecognized: ${unknown.join(", ")}` : ""); return; }

      statusEl.className = "status"; statusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      let geneVals;
      try { geneVals = new Map(); await Promise.all(recs.map(async x => geneVals.set(x.rec.gene_id, await dataset.getGeneValues(x.rec)))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files (CORS?): " + e.message; return; }

      statusEl.textContent = "Computing…";
      const patients = patientsInScope(dataset, cancers);
      let out;
      if (state.mode === "condensed") out = computeCondensed(dims, cancers, recs, geneVals, patients);
      else if (recs.length === 1) out = computeExpandedDims(recs[0], dims, cancers, geneVals, patients);  // 單基因：維度 × 癌種
      else out = computeExpanded(dimById[dimPick.value] || dims[0], cancers, recs, geneVals, patients);    // 多基因：基因 × 癌種

      // 繪圖
      let html = heatmapSVG(out.rows, out.cols, out.cells, { colorMax: out.colorMax, legendLabel: out.legendLabel, caption: out.caption });
      if (out.singleBar) html += out.singleBar;
      resultEl.innerHTML = html;

      // 下載（向量 SVG / 高解析 PNG）
      const dlbar = document.createElement("div"); dlbar.className = "co-toolbar"; dlbar.style.marginTop = "10px";
      dlbar.innerHTML = `<button class="co-mini" id="co-svg">Download SVG (vector)</button><button class="co-mini" id="co-png">Download PNG (high-res)</button>`;
      resultEl.appendChild(dlbar);
      container.querySelector("#co-svg").addEventListener("click", downloadSVG);
      container.querySelector("#co-png").addEventListener("click", downloadPNG);

      const note = document.createElement("div"); note.className = "co-legend";
      note.textContent = `★ FDR q<0.05  ★★<0.01  ★★★<0.001  ★★★★<0.0001  |  grey hatch = no data  |  faded + dot = too few samples, not tested` + (unknown.length ? `  |  Unrecognized: ${unknown.join(", ")}` : "");
      resultEl.appendChild(note);
      lastResults = out.table; $("#co-expres").style.display = "";
      statusEl.textContent = "Done.";
    }

    // ---------- 下載：向量 SVG / 高解析 PNG（投稿用，內容全英文）----------
    function firstSVG() { return resultEl.querySelector("svg"); }
    function triggerDownload(blob, name) { const u = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = u; a.download = name; a.click(); URL.revokeObjectURL(u); }
    function downloadSVG() {
      const svg = firstSVG(); if (!svg) return;
      const xml = new XMLSerializer().serializeToString(svg);
      triggerDownload(new Blob([xml], { type: "image/svg+xml;charset=utf-8" }), `clinical_overview_${state.mode}.svg`);
    }
    function downloadPNG() {
      const svg = firstSVG(); if (!svg) return;
      const vb = svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal : null;
      const W = vb && vb.width ? vb.width : svg.clientWidth, H = vb && vb.height ? vb.height : svg.clientHeight;
      const scale = Math.max(2, Math.ceil(2400 / W));   // 目標寬度約 2400px，足夠投稿
      const xml = new XMLSerializer().serializeToString(svg);
      const img = new Image();
      img.onload = () => {
        const c = document.createElement("canvas"); c.width = Math.round(W * scale); c.height = Math.round(H * scale);
        const ctx = c.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, c.width, c.height);
        ctx.drawImage(img, 0, 0, c.width, c.height);
        c.toBlob(b => triggerDownload(b, `clinical_overview_${state.mode}.png`), "image/png");
      };
      img.onerror = () => alert("PNG export failed (try Download SVG).");
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(xml);
    }

    // ---------- 計算：展開（列=基因、欄=癌種、色=log2FC）----------
    function computeExpanded(D, cancers, recs, geneVals, patients) {
      const a = assignmentOf(D);
      // 每癌種先把病人分桶（只做一次）
      const per = {}; cancers.forEach(c => per[c] = { base: [], adv: [], nField: 0 });
      patients.forEach(p => {
        const cc = per[p.cancer]; if (!cc) return;
        const raw = p.clin[D.field]; if (raw != null && String(raw).trim() !== "") cc.nField++;
        const b = classify(D, raw, a); if (b === "baseline") cc.base.push(p.idx); else if (b === "advanced") cc.adv.push(p.idx);
      });
      const rows = recs.map(r => r.label), cols = cancers;
      const cells = []; const pflat = []; const table = [];
      recs.forEach((r, ri) => {
        const vals = geneVals.get(r.rec.gene_id); const rowCells = [];
        cancers.forEach(c => {
          const cc = per[c]; const base = cc.base.map(i => vals[i]).filter(v => isFinite(v)), adv = cc.adv.map(i => vals[i]).filter(v => isFinite(v));
          const nb = base.length, na = adv.length;
          let cell;
          if (cc.nField === 0) cell = { state: "nodata", tip: `${r.label} / ${c}: no ${D.name} data` };
          else if (nb === 0 || na === 0) cell = { state: "nodata", tip: `${r.label} / ${c}: no contrast (${nb} vs ${na})` };
          else {
            const fc = log2FC(base, adv); const t = mannWhitney(base, adv);
            const hi = Math.max(nb, na), lo = Math.min(nb, na);
            const weak = nb < THR.minPerGroup || na < THR.minPerGroup || hi / lo > THR.maxRatio;
            cell = { value: fc, state: weak ? "weak" : "ok", p: t.p, stars: "", tip: `${r.label} / ${c}: log2FC=${fc.toFixed(2)}, n=${nb}/${na}, p=${t.p.toPrecision(2)}` };
            if (!weak) pflat.push({ ri, ci: cols.indexOf(c), p: t.p });
            table.push({ gene: r.label, cancer: c, n_base: nb, n_adv: na, median_base: median(base), median_adv: median(adv), log2FC: fc, p: t.p });
          }
          rowCells.push(cell);
        });
        cells.push(rowCells);
      });
      // FDR（僅 ok 格）
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = cells[x.ri][cols.indexOf(cols[x.ci])]; cell.q = q[k]; cell.stars = starsOf(q[k]); cell.tip += `, q=${q[k].toPrecision(2)}`; });
      // 把 q 寫回 table
      const qmap = new Map(); pflat.forEach((x, k) => qmap.set(x.ri + "|" + x.ci, q[k]));
      table.forEach(t => { /* q 由 cells tip 已含；table 補 q */ });
      const lbl = labelsOf(D, assignmentOf(D));
      const colorMax = robustMax(cells, 4);
      return {
        rows, cols, cells, colorMax, legendLabel: "log2 FC",
        caption: `${D.name}: ${lbl[1]} vs ${lbl[0]} — log2FC across cancers`,
        table: { mode: "expanded", dim: D.name, baseLabel: lbl[0], advLabel: lbl[1], rows: tableWithQ(table, cells, rows, cols) },
        singleBar: null,
      };
    }

    // ---------- 計算：展開·單基因（列=臨床維度、欄=癌種、色=log2FC）----------
    // 固定一個基因 → 空出一軸，於是 Y=維度、X=癌種，一張圖看該基因在所有維度×癌種的 pattern。
    function computeExpandedDims(geneRec, dims, cancers, geneVals, patients) {
      const vals = geneVals.get(geneRec.rec.gene_id);
      const rows = dims.map(d => d.name), cols = cancers;
      const cells = []; const table = [];
      dims.forEach((D) => {
        const a = assignmentOf(D);
        // 該維度下、各癌種先分桶
        const per = {}; cancers.forEach(c => per[c] = { base: [], adv: [], nField: 0 });
        patients.forEach(p => {
          const cc = per[p.cancer]; if (!cc) return;
          const raw = p.clin[D.field]; if (raw != null && String(raw).trim() !== "") cc.nField++;
          const b = classify(D, raw, a); if (b === "baseline") cc.base.push(p.idx); else if (b === "advanced") cc.adv.push(p.idx);
        });
        const lbl = labelsOf(D, assignmentOf(D));
        const rowCells = []; const pRow = [];
        cancers.forEach((c, ci) => {
          const cc = per[c]; const base = cc.base.map(i => vals[i]).filter(v => isFinite(v)), adv = cc.adv.map(i => vals[i]).filter(v => isFinite(v));
          const nb = base.length, na = adv.length;
          let cell;
          if (cc.nField === 0) cell = { state: "nodata", tip: `${D.name} / ${c}: no data` };
          else if (nb === 0 || na === 0) cell = { state: "nodata", tip: `${D.name} / ${c}: no contrast (${nb} vs ${na})` };
          else {
            const fc = log2FC(base, adv); const t = mannWhitney(base, adv);
            const hi = Math.max(nb, na), lo = Math.min(nb, na);
            const weak = nb < THR.minPerGroup || na < THR.minPerGroup || hi / lo > THR.maxRatio;
            cell = { value: fc, state: weak ? "weak" : "ok", stars: "", tip: `${D.name} / ${c}: ${lbl[1]} vs ${lbl[0]}, log2FC=${fc.toFixed(2)}, n=${nb}/${na}, p=${t.p.toPrecision(2)}` };
            if (!weak) pRow.push({ ci, p: t.p });
            table.push({ gene: geneRec.label, dimension: D.name, cancer: c, baseLabel: lbl[0], advLabel: lbl[1], n_base: nb, n_adv: na, median_base: median(base), median_adv: median(adv), log2FC: fc, p: t.p });
          }
          rowCells.push(cell);
        });
        // 每列（維度）獨立做 FDR
        const q = benjaminiHochberg(pRow.map(x => x.p));
        pRow.forEach((x, k) => { const cell = rowCells[x.ci]; cell.q = q[k]; cell.stars = starsOf(q[k]); cell.tip += `, q=${q[k].toPrecision(2)}`; });
        cells.push(rowCells);
      });
      const colorMax = robustMax(cells, 4);
      return {
        rows, cols, cells, colorMax, legendLabel: "log2 FC",
        caption: `${geneRec.label}: log2FC across cancers (rows = clinical dimensions; advanced vs baseline)`,
        table: { mode: "expanded_dims", rows: tableWithQDims(table, cells, rows, cols) },
        singleBar: null,
      };
    }
    // 單基因展開：q 依「維度列 × 癌種欄」併回 table
    function tableWithQDims(table, cells, rows, cols) {
      return table.map(t => {
        const ri = rows.indexOf(t.dimension), ci = cols.indexOf(t.cancer);
        const q = (cells[ri] && cells[ri][ci] && cells[ri][ci].q != null) ? cells[ri][ci].q : "";
        return { ...t, q };
      });
    }

    // ---------- 計算：濃縮（合併癌種、列=基因、欄=維度兩側、色=z-score）----------
    function computeCondensed(dims, cancers, recs, geneVals, patients) {
      // 欄 = 每個維度兩側
      const cols = []; const colKey = {};
      dims.forEach(D => {
        const lbl = labelsOf(D, assignmentOf(D));
        const cb = { D, side: "baseline", label: lbl[0], idx: [] };
        const ca = { D, side: "advanced", label: lbl[1], idx: [] };
        colKey[D.id + "|baseline"] = cb; colKey[D.id + "|advanced"] = ca;
        cols.push(cb, ca);
      });
      // 分桶（只做一次）
      patients.forEach(p => dims.forEach(D => {
        const b = classify(D, p.clin[D.field], assignmentOf(D));
        if (b) colKey[D.id + "|" + b].idx.push(p.idx);
      }));

      const rows = recs.map(r => r.label);
      const cells = []; const table = [];
      const pPerDimGene = []; // {ri, advCi, p}
      recs.forEach((r, ri) => {
        const vals = geneVals.get(r.rec.gene_id);
        const colVals = cols.map(col => col.idx.map(i => vals[i]).filter(v => isFinite(v)));   // 濾掉缺值樣本（miRNA）
        const colMeans = colVals.map(cv => cv.length ? mean(cv) : null);
        const z = zscoreRow(colMeans);
        const rowCells = cols.map((col, ci) => {
          const nC = colVals[ci].length;
          if (!nC) return { state: "nodata", tip: `${r.label} / ${col.label}: no data` };
          const weak = nC < THR.minPerGroup;
          return { value: z[ci], state: weak ? "weak" : "ok", stars: "", tip: `${r.label} / ${col.label}: z=${(z[ci] ?? 0).toFixed(2)}, n=${nC}, mean=${colMeans[ci].toFixed(2)}` };
        });
        // 每維度成對檢定（star 標在 advanced 欄）
        dims.forEach(D => {
          const cb = colKey[D.id + "|baseline"], ca = colKey[D.id + "|advanced"];
          const advCi = cols.indexOf(ca);
          const bv = cb.idx.map(i => vals[i]).filter(v => isFinite(v)), av = ca.idx.map(i => vals[i]).filter(v => isFinite(v));
          if (bv.length && av.length) {
            const t = mannWhitney(bv, av);
            const weak = bv.length < THR.minPerGroup || av.length < THR.minPerGroup;
            if (!weak) pPerDimGene.push({ ri, advCi, p: t.p });
            table.push({ gene: r.label, dimension: D.name, n_base: bv.length, n_adv: av.length, mean_base: mean(bv), mean_adv: mean(av), p: t.p });
            rowCells[advCi].pairP = t.p;
          }
        });
        cells.push(rowCells);
      });
      const q = benjaminiHochberg(pPerDimGene.map(x => x.p));
      pPerDimGene.forEach((x, k) => { const cell = cells[x.ri][x.advCi]; if (cell.state !== "nodata") { cell.q = q[k]; cell.stars = starsOf(q[k]); cell.tip += `, q=${q[k].toPrecision(2)}`; } });

      const colorMax = robustMax(cells, 3);
      const out = {
        rows, cols: cols.map(c => c.label), cells, colorMax, legendLabel: "row z-score",
        caption: `Row z-score (merged ${cancers.length} cancers)`,
        table: { mode: "condensed", rows: table },
        singleBar: null,
      };

      // 單基因：附 image 2 那種 bar（各組原始 TPM mean±SD）
      if (recs.length === 1) {
        const vals = geneVals.get(recs[0].rec.gene_id);
        const unit = recs[0].rec.assay === "mirna_rpm" ? "RPM" : "TPM";   // miRNA 標 RPM
        const groups = cols.map((col, ci) => {
          const gv = col.idx.map(i => vals[i]).filter(v => isFinite(v));
          if (!gv.length) return { name: col.label, state: "nodata" };
          const weak = gv.length < THR.minPerGroup;
          const cell = cells[0][ci];
          const nG = gv.length, sem = nG > 1 ? sd(gv) / Math.sqrt(nG) : 0;
          return { name: col.label, m: mean(gv), err: sem, color: col.side === "baseline" ? "#cbd5e1" : "#64748b", state: weak ? "weak" : "ok", stars: cell.stars };
        });
        out.singleBar = multiBarSVG(groups, { caption: `${recs[0].label} expression (mean ± SEM)`, ylabel: unit });
      }
      return out;
    }

    function robustMax(cells, cap) {
      const vals = [];
      cells.forEach(row => row.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) vals.push(Math.abs(c.value)); }));
      if (!vals.length) return 1;
      const m = Math.max(...vals);
      return Math.min(cap, Math.max(0.5, Math.ceil(m * 10) / 10));
    }

    // 把 q 併回展開模式的 table（依 cells 的 q）
    function tableWithQ(table, cells, rows, cols) {
      return table.map(t => {
        const ri = rows.indexOf(t.gene), ci = cols.indexOf(t.cancer);
        const q = (cells[ri] && cells[ri][ci] && cells[ri][ci].q != null) ? cells[ri][ci].q : "";
        return { ...t, q };
      });
    }

    // ---------- 匯出結果表 ----------
    function exportResults() {
      if (!lastResults) return;
      const rows = lastResults.rows; if (!rows.length) return;
      const cols = Object.keys(rows[0]);
      const esc = v => (v == null ? "" : (typeof v === "number" ? (Number.isInteger(v) ? v : v.toFixed(4)) : String(v)));
      let csv = cols.join(",") + "\n";
      rows.forEach(r => { csv += cols.map(c => esc(r[c])).join(",") + "\n"; });
      const blob = new Blob([csv], { type: "text/csv" }); const url = URL.createObjectURL(blob);
      const a = document.createElement("a"); a.href = url; a.download = `clinical_overview_${state.mode}.csv`; a.click(); URL.revokeObjectURL(url);
    }

    // GOIs 與其他分析（Survival）共享：別處改動時同步本頁
    onGOIsChanged(list => { state.genes = list.slice(); if (document.activeElement !== genesEl) genesEl.value = list.join(", "); commit(); refreshDimPick(); });

    // ---------- 初始化 ----------
    syncGenes(); refreshNamed(); renderAll(); refreshDimPick();
  },
};
