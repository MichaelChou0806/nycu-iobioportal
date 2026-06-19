// =====================================================================
// analyses/survivalGroups.js —「Survival · Custom Groups」（OS）
// ---------------------------------------------------------------------
// 假設驅動的組合分組生存分析（有別於基礎 Survival 的掃描式）。
//   使用者用「幾個基因的 high/low 組合」自訂兩組病人（A / B），再比生存。
//   例：HH vs LL、HH vs rest、LL vs rest，或逐格自訂。
//
// 第一版範圍（刻意聚焦，能跑能驗證）：
//   2 基因 × 單癌別 × median 二分 → 一條 KM（A vs B）+ log-rank p + Cox HR
//   每組 n、每格 signature 的 n 都顯示。
//
// 架構預留（升級不重寫，只動 UI）：
//   · 資料結構以 N 基因通用：sigKeys(n) 產生 2ⁿ 個組合鍵，assignment 用字典。
//   · 二分抽象為 splitGroups(method)，與基礎 Survival 共用；先鎖 median。
//   · 病人取得後先 filter 的位置已留好 → 未來插臨床 subset。
//   · 透過 Dataset 介面取資料 → dataset-agnostic，OSCC 可直接套。
//   未來：臨床 subset、multi-cancer（per-cancer / pooled）、cutoff 選擇、多基因 UI。
//
// 統計核心（kaplanMeier / logRank / coxPH1）見 stats.js（已驗證）。
// =====================================================================

import { patientsInScope, classify, loadDimensions } from "../core/dimensions.js";
import { kaplanMeier, logRank, coxPH1, logRankStratified, coxPH1Stratified, benjaminiHochberg, pStars } from "../core/stats.js";
import { kmCurveSVG, heatmapSVG } from "../core/plots.js";

const DAYS_PER_MONTH = 30.4375;
const A_COLOR = "#ef4444", B_COLOR = "#3b82f6";   // A=紅, B=藍

let stylesInjected = false;
function injectStyles() {
  if (stylesInjected) return; stylesInjected = true;
  const css = `
    .sg-sec{margin-bottom:18px}
    .sg-intro{font-size:12.5px;color:var(--muted);line-height:1.5;margin-bottom:12px}
    .sg-h3{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;margin:0 0 8px}
    .sg-sum{font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em;font-weight:600;cursor:pointer;padding:4px 0}
    .sg-toolbar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px;align-items:center}
    .sg-ctl{font-size:11.5px;color:var(--muted)}
    .sg-ctl input{padding:6px 9px;font-size:13px;margin-left:4px;border:1px solid var(--line);border-radius:6px;width:130px}
    .sg-mini{padding:4px 10px;font-size:12px;border-radius:6px;background:#fff;color:var(--accent);border:1px solid var(--accent);cursor:pointer}
    .sg-run-big{padding:9px 20px;font-weight:600}
    .sg-sep{width:1px;height:18px;background:var(--line);margin:0 4px}
    .sg-chips{display:flex;flex-wrap:wrap;gap:6px}
    .sg-chip{padding:5px 11px;border:1px solid var(--line);border-radius:20px;background:#fff;font-size:12px;cursor:grab;user-select:none}
    .sg-chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .sg-chip.dragging{opacity:.4}
    .sg-presets{align-items:center}
    .sg-presets-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em;margin-right:2px}
    .sg-grid{display:flex;flex-wrap:wrap;gap:10px}
    .sg-cell{border:1px solid var(--line);border-radius:8px;padding:10px 12px;min-width:120px}
    .sg-sig{font:600 15px monospace;letter-spacing:.06em;margin-bottom:6px;color:#1f2733}
    .sg-seg{display:inline-flex;border:1px solid var(--line);border-radius:6px;overflow:hidden}
    .sg-seg button{padding:4px 11px;font-size:12px;border:none;background:#fff;cursor:pointer;border-left:1px solid var(--line)}
    .sg-seg button:first-child{border-left:none}
    .sg-seg button.onA{background:#ef4444;color:#fff}
    .sg-seg button.onB{background:#3b82f6;color:#fff}
    .sg-seg button.onX{background:#9aa3af;color:#fff}
    .sg-n{font-size:11px;color:var(--muted);margin-top:6px;min-height:14px}
    .sg-note{font-size:11.5px;color:var(--muted)}
    .sg-result{overflow-x:auto}
    .sg-result svg{max-width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:8px;margin-top:10px}
    .sg-legend{font-size:11.5px;color:var(--muted);margin-top:6px;line-height:1.5}
    .sg-divider{height:1px;background:var(--line);margin:22px 0}
    .sg-sub{font-weight:400;text-transform:none;letter-spacing:0;font-size:11.5px;color:var(--muted)}
    .sg-scr-intro{font-size:12px;color:var(--muted);line-height:1.5;margin-bottom:10px}
    .sg-screen-dims{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:12px}
    .sg-dim{display:flex;align-items:center;gap:6px;padding:6px 11px;border:1px solid var(--line);border-radius:6px;font-size:12.5px;cursor:pointer;user-select:none}
    .sg-dim.on{background:var(--accent);color:#fff;border-color:var(--accent)}
    .sg-dim input{cursor:pointer;margin:0}
    .sg-stale{font-size:11.5px;color:#b45309;padding:4px 0}
  `;
  const el = document.createElement("style"); el.textContent = css; document.head.appendChild(el);
}

const LS = "tcga-tool:survivalGroups";
function loadState() { try { return JSON.parse(localStorage.getItem(LS)) || {}; } catch (e) { return {}; } }
function saveState(st) { try { localStorage.setItem(LS, JSON.stringify(st)); } catch (e) {} }
function moveInArray(arr, fromId, toId) {
  if (fromId == null || fromId === toId) return arr.slice();
  const a = arr.slice(); const fi = a.indexOf(fromId); if (fi < 0) return a;
  a.splice(fi, 1); let ti = a.indexOf(toId); if (ti < 0) ti = a.length; a.splice(ti, 0, fromId); return a;
}
function quantileSorted(sorted, f) { const p = f * (sorted.length - 1), lo = Math.floor(p), hi = Math.ceil(p); return sorted[lo] + (sorted[hi] - sorted[lo]) * (p - lo); }
// 切法分組：1=High, 0=Low, -1=丟棄（與基礎 Survival 同一套；本版先鎖 median）
function splitGroups(vals, method) {
  const sorted = [...vals].sort((a, b) => a - b);
  const frac = method === "tertile" ? 1 / 3 : method === "quartile" ? 1 / 4 : 0.5;
  const lowCut = quantileSorted(sorted, frac), highCut = quantileSorted(sorted, 1 - frac);
  return vals.map(v => method === "median" ? (v > highCut ? 1 : 0) : (v <= lowCut ? 0 : (v > highCut ? 1 : -1)));
}
function fmtP(p) { return p < 0.001 ? "< 0.001" : p.toFixed(3); }

// N 基因通用：產生所有 2ⁿ 個 signature 鍵，如 n=2 → ["HH","HL","LH","LL"]
function sigKeys(n) { let out = [""]; for (let i = 0; i < n; i++) out = out.flatMap(s => [s + "H", s + "L"]); return out; }

// 維度兩側的人類可讀標籤（複用 clinicalOverview 邏輯）
function labelsOf(dim, assignment) {
  if (dim.type === "binary") {
    if (dim.numericSplit) { const c = (assignment && assignment.cutoff != null) ? assignment.cutoff : dim.numericSplit.cutoff; return [`<${c}`, `>=${c}`]; }
    if (dim.baseline && dim.advanced) return [dim.baseline.label, dim.advanced.label];
  }
  if (dim.type === "ordinal") return [ordinalRange(dim, assignment, "baseline"), ordinalRange(dim, assignment, "advanced")];
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
// 用維度定義的 default 產生 classify 用的 assignment（subset 不讓使用者覆寫級別）
function assignmentFor(dim) {
  if (dim.type === "ordinal") { const a = {}; dim.levels.forEach(lv => a[lv.id] = lv.default); return a; }
  if (dim.numericSplit) return { cutoff: dim.numericSplit.cutoff };
  return null;
}
function splitLabel(m) { return m === "tertile" ? "tertile (top/bottom 1/3)" : m === "quartile" ? "quartile (Q1 vs Q4)" : "median split"; }

export const survivalGroups = {
  id: "survivalGroups",
  name: "Advanced Survival",

  async mount(container, { dataset }) {
    injectStyles();
    // 載入臨床維度定義（subset 用；失敗則此功能降級，不影響核心 KM）
    let dims = [];
    try { const dimUrl = (dataset.config && dataset.config.dimensionsUrl) || "config/dimensions.tcga.json"; const DEF = await loadDimensions(dimUrl); dims = (DEF && DEF.dimensions) || []; } catch (e) { dims = []; }
    const dimById = Object.fromEntries(dims.map(d => [d.id, d]));
    const availCancers = dataset.cancers.map(c => c.code);
    const nTumorOf = Object.fromEntries(dataset.cancers.map(c => [c.code, c.n_tumor]));
    const saved = loadState();
    const state = {
      nGenes: (saved.nGenes >= 1 && saved.nGenes <= 4) ? saved.nGenes : 2,
      genes: Array.isArray(saved.genes) ? saved.genes.slice(0, 4) : [],
      cancers: (Array.isArray(saved.cancers) ? saved.cancers.filter(c => availCancers.includes(c)) : []).concat(availCancers.filter(c => !(saved.cancers || []).includes(c))),
      selectedCancers: Array.isArray(saved.selectedCancers) ? saved.selectedCancers.filter(c => availCancers.includes(c)) : [],
      split: ["median", "tertile", "quartile"].includes(saved.split) ? saved.split : "median",
      months: saved.months != null ? saved.months : 0,
      // signature → "A" | "B" | "exclude"；預設 HH vs LL
      assign: saved.assign && typeof saved.assign === "object" ? saved.assign : { HH: "A", LL: "B", HL: "exclude", LH: "exclude" },
      subset: saved.subset || "",         // ""=無 subset；否則 "dimId|baseline" / "dimId|advanced"
      view: saved.view === "pooled" ? "pooled" : "percancer",   // 多癌別：percancer | pooled
      scheme: saved.scheme === "rb" ? "rb" : "rg",
      colorMax: saved.colorMax != null ? saved.colorMax : 0,    // 0=auto
      // ── Subset screening（各臨床亞組各自 A vs B 的 2D heatmap）──
      selectedScreenDims: Array.isArray(saved.selectedScreenDims) ? saved.selectedScreenDims : [],
      scrScheme: saved.scrScheme === "rb" ? "rb" : "rg",
      scrColorMax: saved.scrColorMax != null ? saved.scrColorMax : 0,
    };
    const subsetOptions = `<option value="">No clinical subset (all tumors)</option>` +
      dims.map(d => { const [bl, al] = labelsOf(d, assignmentFor(d)); const nm = d.name || d.id; return `<option value="${d.id}|baseline">${nm}: ${bl}</option><option value="${d.id}|advanced">${nm}: ${al}</option>`; }).join("");
    let dragCancer = null, lastSVGName = "km_groups", lastKM = null, scrLastSVGName = "screening";

    container.innerHTML = `
      <div class="card">
        <div class="sg-sec">
          <div class="sg-intro">Define custom patient groups by combining two genes' high/low levels, then compare survival (OS). Optionally restrict to a clinical subset first (e.g. node-positive only). Pick which signatures go to group A, group B, or are excluded. Each group's n is shown. This is the hypothesis-driven counterpart to the scanning-style Survival tab.</div>
        </div>

        <div class="sg-sec">
          <h3 class="sg-h3">Genes</h3>
          <div class="sg-toolbar">
            <label class="sg-ctl"># of genes<select id="sg-ngenes">
              <option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option>
            </select></label>
            <span class="sg-ctl" style="opacity:.7">→ 2ⁿ signature groups below</span>
          </div>
          <div id="sg-genes" class="sg-toolbar" style="margin-top:8px"></div>
        </div>

        <details class="sg-sec" open>
          <summary class="sg-sum">Cancers — pick one for a KM curve, several for heatmap / pooled</summary>
          <div class="sg-toolbar" style="margin-top:8px">
            <button class="sg-mini" id="sg-alpha">Sort A–Z</button>
            <button class="sg-mini" id="sg-byn">Sort by N</button>
            <button class="sg-mini" id="sg-def">Default order</button>
            <span class="sg-sep"></span>
            <button class="sg-mini" id="sg-all">Select all</button>
            <button class="sg-mini" id="sg-none">Clear</button>
          </div>
          <div id="sg-cancers" class="sg-chips"></div>
        </details>

        <div class="sg-sec">
          <h3 class="sg-h3">Groups</h3>
          <div class="sg-toolbar">
            <label class="sg-ctl">Split<select id="sg-split">
              <option value="median">Median (50/50)</option>
              <option value="tertile">Tertile (top/bottom 1/3)</option>
              <option value="quartile">Quartile (Q4 vs Q1)</option>
            </select></label>
            <label class="sg-ctl">Follow-up months<input type="number" id="sg-months" min="0" step="1" style="width:72px" placeholder="all"></label>
          </div>
          <div class="sg-toolbar sg-presets">
            <span class="sg-presets-label">Presets</span>
            <button class="sg-mini" id="sg-p-hhll">All-High vs All-Low</button>
            <button class="sg-mini" id="sg-p-hhrest">All-High vs rest</button>
            <button class="sg-mini" id="sg-p-llrest">All-Low vs rest</button>
          </div>
          <div id="sg-grid" class="sg-grid"></div>
          <div class="sg-note" style="margin-top:8px">Tumor only · endpoint OS. Signature reads genes in the order entered (with 2 genes, HL = Gene 1 High, Gene 2 Low). With 1 gene it's simply High vs Low.</div>
        </div>

        <div class="sg-divider"></div>

        <div class="sg-sec">
          <h3 class="sg-h3">Single comparison <span class="sg-sub">— one KM (or per-cancer HR), optionally within one clinical subset</span></h3>
          <div class="sg-toolbar">
            <label class="sg-ctl">Clinical subset<select id="sg-subset" style="min-width:220px">${subsetOptions}</select></label>
            <label class="sg-ctl" id="sg-view-wrap">Multi-cancer<select id="sg-view">
              <option value="percancer">Per-cancer (HR heatmap)</option>
              <option value="pooled">Pooled (one KM, stratified)</option>
            </select></label>
            <label class="sg-ctl" id="sg-scheme-wrap">Colors<select id="sg-scheme">
              <option value="rg">Red–Green</option><option value="rb">Red–Blue</option>
            </select></label>
            <label class="sg-ctl" id="sg-cmax-wrap">Scale max<input type="number" id="sg-cmax" min="0" step="0.1" style="width:64px" placeholder="auto"></label>
          </div>
          <div class="sg-toolbar">
            <button id="sg-run">Run comparison</button>
            <span class="sg-sep"></span>
            <button class="sg-mini" id="sg-prism" style="display:none">Export CSV</button>
            <button class="sg-mini" id="sg-svg" style="display:none">SVG</button>
            <button class="sg-mini" id="sg-png" style="display:none">PNG</button>
          </div>
          <div id="sg-status" class="status"></div>
          <div id="sg-result" class="sg-result"></div>
        </div>

        <div class="sg-divider"></div>

        <div class="sg-sec">
          <h3 class="sg-h3">Subset screening <span class="sg-sub">— scan the A-vs-B effect across clinical subgroups (exploratory)</span></h3>
          <div class="sg-scr-intro">Pick clinical variables to screen. Each variable's two sides become columns (e.g. N stage → N0 and N+); rows are your selected cancers plus a Pooled row. Each cell is the A-vs-B log2 HR within that subgroup. Exploratory — it suggests where the effect is stronger, but does not by itself test for interaction.</div>
          <div class="sg-toolbar sg-presets">
            <span class="sg-presets-label">Variables</span>
            <button class="sg-mini" id="sg-scr-all">Select all</button>
            <button class="sg-mini" id="sg-scr-clear">Clear</button>
          </div>
          <div id="sg-screen-dims" class="sg-screen-dims"></div>
          <div class="sg-toolbar">
            <label class="sg-ctl">Colors<select id="sg-scr-scheme">
              <option value="rg">Red–Green</option><option value="rb">Red–Blue</option>
            </select></label>
            <label class="sg-ctl">Scale max<input type="number" id="sg-scr-cmax" min="0" step="0.1" style="width:64px" placeholder="auto"></label>
          </div>
          <div class="sg-toolbar">
            <button id="sg-run-screen" class="sg-run-big">Run screening</button>
            <span class="sg-sep"></span>
            <button class="sg-mini" id="sg-scr-svg" style="display:none">SVG</button>
            <button class="sg-mini" id="sg-scr-png" style="display:none">PNG</button>
          </div>
          <div id="sg-scr-status" class="status"></div>
          <div id="sg-scr-result" class="sg-result"></div>
        </div>
      </div>`;

    const $ = s => container.querySelector(s);
    const monthsEl = $("#sg-months"), splitSel = $("#sg-split"), subsetSel = $("#sg-subset"), ngenesSel = $("#sg-ngenes"), genesBox = $("#sg-genes"), viewSel = $("#sg-view"), schemeSel = $("#sg-scheme"), cmaxEl = $("#sg-cmax");
    const chipBox = $("#sg-cancers"), gridBox = $("#sg-grid"), statusEl = $("#sg-status"), resultEl = $("#sg-result");
    const scrSchemeSel = $("#sg-scr-scheme"), scrCmaxEl = $("#sg-scr-cmax"), screenDimsBox = $("#sg-screen-dims"), scrStatusEl = $("#sg-scr-status"), scrResultEl = $("#sg-scr-result");
    function commit() { saveState(state); }

    if (state.months > 0) monthsEl.value = state.months;
    splitSel.value = state.split; subsetSel.value = state.subset; ngenesSel.value = state.nGenes;
    viewSel.value = state.view; schemeSel.value = state.scheme; if (state.colorMax > 0) cmaxEl.value = state.colorMax;
    monthsEl.addEventListener("change", () => { const v = Number(monthsEl.value); state.months = (monthsEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); markStale(); });
    splitSel.addEventListener("change", () => { state.split = splitSel.value; commit(); markStale(); });
    subsetSel.addEventListener("change", () => { state.subset = subsetSel.value; commit(); });
    viewSel.addEventListener("change", () => { state.view = viewSel.value; commit(); updateControls(); });
    schemeSel.addEventListener("change", () => { state.scheme = schemeSel.value; commit(); });
    cmaxEl.addEventListener("change", () => { const v = Number(cmaxEl.value); state.colorMax = (cmaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); });
    // ── Subset screening 控制項 ──
    scrSchemeSel.value = state.scrScheme; if (state.scrColorMax > 0) scrCmaxEl.value = state.scrColorMax;
    scrSchemeSel.addEventListener("change", () => { state.scrScheme = scrSchemeSel.value; commit(); });
    scrCmaxEl.addEventListener("change", () => { const v = Number(scrCmaxEl.value); state.scrColorMax = (scrCmaxEl.value === "" || !isFinite(v) || v <= 0) ? 0 : v; commit(); });
    // 共享設定改變 → 兩邊結果都失效（避免新舊並存誤導）
    function markStale() {
      if (resultEl.querySelector("svg") || scrResultEl.querySelector("svg")) {
        resultEl.innerHTML = ""; scrResultEl.innerHTML = "";
        statusEl.className = "status sg-stale"; statusEl.textContent = "Shared settings changed — re-run.";
        scrStatusEl.className = "status sg-stale"; scrStatusEl.textContent = "Shared settings changed — re-run.";
        ["sg-prism", "sg-svg", "sg-png", "sg-scr-svg", "sg-scr-png"].forEach(id => $("#" + id).style.display = "none");
      }
    }

    // 動態渲染基因輸入框（數量 = state.nGenes）
    function renderGenes() {
      genesBox.innerHTML = "";
      for (let i = 0; i < state.nGenes; i++) {
        const lab = document.createElement("label"); lab.className = "sg-ctl"; lab.textContent = `Gene ${i + 1}`;
        const inp = document.createElement("input"); inp.value = state.genes[i] || ""; inp.placeholder = i === 0 ? "e.g. KDELR1" : i === 1 ? "e.g. KDELR2" : "gene symbol";
        inp.addEventListener("input", () => { state.genes[i] = inp.value.trim(); commit(); markStale(); });
        lab.appendChild(inp); genesBox.appendChild(lab);
      }
    }
    // 基因數變 → 截斷基因清單、把分配重置為 All-High vs All-Low、重畫輸入框與格子
    ngenesSel.addEventListener("change", () => {
      state.nGenes = Number(ngenesSel.value);
      state.genes = state.genes.slice(0, state.nGenes);
      applyPreset("hhll");        // 內含 renderGrid + commit + markStale
      renderGenes();
    });

    // ---- 癌別 chips（複用 Survival 的互動；本版單癌別）----
    function renderCancers() {
      chipBox.innerHTML = "";
      state.cancers.forEach(code => {
        const chip = document.createElement("div");
        chip.className = "sg-chip" + (state.selectedCancers.includes(code) ? " on" : "");
        chip.dataset.code = code; chip.textContent = code; chip.draggable = true;
        chip.addEventListener("click", () => {
          if (state.selectedCancers.includes(code)) state.selectedCancers = state.selectedCancers.filter(c => c !== code);
          else state.selectedCancers.push(code);
          commit(); renderCancers(); updateControls(); markStale();
        });
        chip.addEventListener("dragstart", () => { dragCancer = code; chip.classList.add("dragging"); });
        chip.addEventListener("dragend", () => chip.classList.remove("dragging"));
        chip.addEventListener("dragover", e => e.preventDefault());
        chip.addEventListener("drop", e => { e.preventDefault(); state.cancers = moveInArray(state.cancers, dragCancer, code); commit(); renderCancers(); });
        chipBox.appendChild(chip);
      });
    }
    $("#sg-alpha").addEventListener("click", () => { state.cancers = [...state.cancers].sort(); commit(); renderCancers(); });
    $("#sg-byn").addEventListener("click", () => { state.cancers = [...state.cancers].sort((a, b) => (nTumorOf[b] || 0) - (nTumorOf[a] || 0)); commit(); renderCancers(); });
    $("#sg-def").addEventListener("click", () => { state.cancers = availCancers.slice(); commit(); renderCancers(); });
    $("#sg-all").addEventListener("click", () => { state.selectedCancers = state.cancers.slice(); commit(); renderCancers(); updateControls(); markStale(); });
    $("#sg-none").addEventListener("click", () => { state.selectedCancers = []; commit(); renderCancers(); updateControls(); markStale(); });

    // ---- 分組分配格（N 基因通用，本版渲染 2 基因 = 4 格）----
    function ensureAssignKeys() {                       // 補齊/清掉鍵，讓 assign 對齊目前基因數
      const keys = sigKeys(state.nGenes);
      const next = {};
      keys.forEach(k => next[k] = state.assign[k] || "exclude");
      state.assign = next;
    }
    function renderGrid(counts) {                       // counts: {sig: n}（Run 後才有）
      ensureAssignKeys();
      const keys = sigKeys(state.nGenes);
      gridBox.innerHTML = keys.map(k => {
        const a = state.assign[k];
        const nTxt = counts && counts[k] != null ? `n = ${counts[k]}` : "";
        return `<div class="sg-cell">
          <div class="sg-sig">${k}</div>
          <div class="sg-seg" data-key="${k}">
            <button data-v="A" class="${a === 'A' ? 'onA' : ''}">A</button>
            <button data-v="B" class="${a === 'B' ? 'onB' : ''}">B</button>
            <button data-v="exclude" class="${a === 'exclude' ? 'onX' : ''}">—</button>
          </div>
          <div class="sg-n">${nTxt}</div>
        </div>`;
      }).join("");
      gridBox.querySelectorAll(".sg-seg button").forEach(btn => {
        btn.addEventListener("click", () => {
          const key = btn.parentElement.dataset.key, v = btn.dataset.v;
          state.assign[key] = v; commit(); renderGrid(counts); markStale();
        });
      });
    }
    function applyPreset(name) {
      const keys = sigKeys(state.nGenes);
      const allH = "H".repeat(state.nGenes), allL = "L".repeat(state.nGenes);
      const next = {};
      if (name === "hhll") keys.forEach(k => next[k] = k === allH ? "A" : k === allL ? "B" : "exclude");
      else if (name === "hhrest") keys.forEach(k => next[k] = k === allH ? "A" : "B");
      else if (name === "llrest") keys.forEach(k => next[k] = k === allL ? "A" : "B");
      state.assign = next; commit(); renderGrid(); markStale();
    }
    $("#sg-p-hhll").addEventListener("click", () => applyPreset("hhll"));
    $("#sg-p-hhrest").addEventListener("click", () => applyPreset("hhrest"));
    $("#sg-p-llrest").addEventListener("click", () => applyPreset("llrest"));

    // 該組包含哪些 signature → 顯示用標籤
    function groupLabel(which) {
      const all = sigKeys(state.nGenes);
      const thisK = all.filter(k => state.assign[k] === which);
      if (!thisK.length) return "(none)";
      if (thisK.length === 1) return thisK[0];                          // 單一 signature → 直接顯示（HH）
      const otherW = which === "A" ? "B" : "A";
      const otherK = all.filter(k => state.assign[k] === otherW);
      if (thisK.length + otherK.length === all.length) return "others"; // 是另一組的補集 → others
      return thisK.join(" / ");                                         // 否則列出（HL / LH）
    }

    // ---- 共用：subset filter / cohort / 分組 / KM 渲染（單癌別、per-cancer、pooled 都用）----
    function applySubset(patients) {
      if (!state.subset) return { patients, txt: "" };
      const [dimId, bucket] = state.subset.split("|");
      const dim = dimById[dimId];
      if (!dim) return { patients, txt: "" };
      const asg = assignmentFor(dim); const [bl, al] = labelsOf(dim, asg);
      return { patients: patients.filter(p => classify(dim, p.clin[dim.field], asg) === bucket), txt: `${dim.name || dim.id} = ${bucket === "baseline" ? bl : al}` };
    }
    // 病人 → cohort（只留有效 expr + 生存；套月數截斷）
    function buildCohort(patients, geneVals) {
      const cohort = [];
      for (const p of patients) {
        const os = Number(p.clin["OS"]), t = Number(p.clin["OS.time"]);
        if (!isFinite(os) || !isFinite(t) || t < 0) continue;
        const exprs = geneVals.map(v => v[p.idx]);
        if (exprs.some(x => !isFinite(x))) continue;
        let tm = t / DAYS_PER_MONTH, e = os ? 1 : 0;
        if (state.months > 0 && tm > state.months) { e = 0; tm = state.months; }
        cohort.push({ exprs, tm, e });
      }
      return cohort;
    }
    // cohort → {A, B, counts}（每基因二分 + signature + 分配）；二分在「傳入的 cohort 內」算分位
    function cohortToAB(cohort, nGene) {
      const levels = []; for (let gi = 0; gi < nGene; gi++) levels.push(splitGroups(cohort.map(c => c.exprs[gi]), state.split));
      const counts = {}; sigKeys(nGene).forEach(k => counts[k] = 0);
      const A = [], B = [];
      cohort.forEach((c, k) => {
        const bits = []; for (let gi = 0; gi < nGene; gi++) bits.push(levels[gi][k]);
        if (bits.includes(-1)) return;                         // 非 median 時被切到中間 → 丟棄
        const key = bits.map(b => b === 1 ? "H" : "L").join("");
        counts[key]++;
        const grp = state.assign[key];
        if (grp === "A") A.push(c); else if (grp === "B") B.push(c);
      });
      return { A, B, counts };
    }
    // A vs B → KM（strata 有值用癌別分層統計）；回傳是否成功
    function renderKM_AB(A, B, strata, caption, xCap) {
      if (!A.length || !B.length) return false;
      const kmA = kaplanMeier(A.map(r => r.tm), A.map(r => r.e));
      const kmB = kaplanMeier(B.map(r => r.tm), B.map(r => r.e));
      const allTm = A.concat(B).map(r => r.tm), allE = A.concat(B).map(r => r.e), allG = A.map(() => 1).concat(B.map(() => 0));
      let lr, cox, strat = false;
      if (strata) { lr = logRankStratified(allTm, allE, allG, strata); cox = coxPH1Stratified(allTm, allE, allG, strata); strat = true; }
      else { lr = logRank(allTm, allE, allG); cox = coxPH1(allTm, allE, allG); }
      const curves = [{ label: groupLabel("A"), color: A_COLOR, km: kmA, n: A.length }, { label: groupLabel("B"), color: B_COLOR, km: kmB, n: B.length }];
      const hrText = isFinite(cox.hr) ? `HR (A vs B${strat ? ", stratified" : ""}) = ${cox.hr.toFixed(2)} (95% CI ${cox.ciLow.toFixed(2)}–${cox.ciHigh.toFixed(2)})` : "HR = n/a";
      resultEl.innerHTML = kmCurveSVG(curves, { caption, xLabel: xCap, pText: `Log-rank p = ${fmtP(lr.p)}${strat ? " (stratified)" : ""}`, hrText, xMax: state.months > 0 ? state.months : undefined });
      return true;
    }
    // 控制項顯示：多癌別才有 view；配色/scale 只在 per-cancer heatmap
    function show(sel, on) { const el = $(sel); if (el) el.style.display = on ? "" : "none"; }
    function updateControls() {
      const nC = state.selectedCancers.length;
      const isHeatmap = nC > 1 && state.view === "percancer";
      show("#sg-view-wrap", nC > 1);
      show("#sg-scheme-wrap", isHeatmap);
      show("#sg-cmax-wrap", isHeatmap);
    }

    // ---- RUN ----
    async function run() {
      resultEl.innerHTML = "";
      ["sg-prism", "sg-svg", "sg-png"].forEach(id => $("#" + id).style.display = "none");
      // 1. 解析兩個基因
      const raw = state.genes.slice(0, state.nGenes);
      if (raw.length < state.nGenes || raw.some(g => !g)) { statusEl.className = "status err"; statusEl.textContent = `Enter ${state.nGenes} gene(s).`; return; }
      const recs = []; const bad = [];
      raw.forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) bad.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (bad.length) { statusEl.className = "status err"; statusEl.textContent = "Unrecognized gene(s): " + bad.join(", "); return; }
      // 2. 至少要選一個癌別
      if (!state.selectedCancers.length) { statusEl.className = "status err"; statusEl.textContent = "Select at least one cancer."; return; }
      const cancers = state.cancers.filter(c => state.selectedCancers.includes(c));   // 依排序保序
      // 3. 至少 A、B 各要有 signature
      const aKeys = sigKeys(state.nGenes).filter(k => state.assign[k] === "A");
      const bKeys = sigKeys(state.nGenes).filter(k => state.assign[k] === "B");
      if (!aKeys.length || !bKeys.length) { statusEl.className = "status err"; statusEl.textContent = "Assign at least one signature to A and one to B."; return; }

      statusEl.className = "status"; statusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      let geneVals;
      try { geneVals = await Promise.all(recs.map(x => dataset.getGeneValues(x.rec))); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "Failed to load gene files (CORS?): " + e.message; return; }

      statusEl.textContent = "Computing…";
      const gLabel = recs.map(r => r.label).join("·");
      const xCap = state.months > 0 ? `Months (capped ${state.months})` : "Months";
      let isKM = false, ok = true;

      if (cancers.length === 1) {
        // ---- 單癌別：一條 KM(A vs B) ----
        const cancer = cancers[0];
        const { patients, txt } = applySubset(patientsInScope(dataset, [cancer]));
        const cohort = buildCohort(patients, geneVals);
        if (cohort.length < 4) { statusEl.className = "status err"; statusEl.textContent = `Too few patients with valid data (n=${cohort.length}).`; return; }
        const { A, B, counts } = cohortToAB(cohort, state.nGenes);
        renderGrid(counts);
        if (!A.length || !B.length) { statusEl.className = "status err"; statusEl.textContent = `One group is empty (A n=${A.length}, B n=${B.length}). Adjust assignment.`; return; }
        lastSVGName = `KM_${recs.map(r => r.label).join("_")}_${cancer}`;
        renderKM_AB(A, B, null, `${gLabel}\n${cancer}${state.months > 0 ? ` · ${state.months}-mo OS` : ""}${txt ? " · " + txt : ""}`, xCap);
        lastKM = { A, B }; isKM = true;
        const evA = A.reduce((s, x) => s + x.e, 0), evB = B.reduce((s, x) => s + x.e, 0);
        const discarded = cohort.length - Object.values(counts).reduce((a, b) => a + b, 0);
        const note = document.createElement("div"); note.className = "sg-legend";
        note.textContent = `${groupLabel("A")}: n=${A.length}, events=${evA} · ${groupLabel("B")}: n=${B.length}, events=${evB}${txt ? ` · subset: ${txt}` : ""} · OS · ${splitLabel(state.split)}${discarded > 0 ? ` · middle dropped: n=${discarded}` : ""}${state.months > 0 ? ` · capped ${state.months} mo` : ""}`;
        resultEl.appendChild(note);
      } else if (state.view === "pooled") {
        // ---- 多癌別 pooled：各癌別內分組、癌別分層的一條 KM ----
        renderGrid();
        ok = drawPooledKM(cancers, geneVals, recs, gLabel, xCap);
        isKM = ok;
      } else {
        // ---- 多癌別 per-cancer：每癌別 A vs B 的 HR heatmap ----
        renderGrid();
        drawHeatmap(cancers, geneVals, recs);
      }
      if (!ok) return;

      if (isKM) $("#sg-prism").style.display = "";
      $("#sg-svg").style.display = ""; $("#sg-png").style.display = "";
      statusEl.textContent = "Done.";
    }

    // ---- 多癌別 pooled：每癌別內各自二分→分組，合併後按癌別分層 ----
    function drawPooledKM(cancers, geneVals, recs, gLabel, xCap) {
      const A = [], B = [], aStrata = [], bStrata = []; let subsetTxt = "";
      cancers.forEach(c => {
        const { patients, txt } = applySubset(patientsInScope(dataset, [c])); subsetTxt = txt;
        const { A: a, B: b } = cohortToAB(buildCohort(patients, geneVals), state.nGenes);
        a.forEach(r => { A.push(r); aStrata.push(c); });
        b.forEach(r => { B.push(r); bStrata.push(c); });
      });
      if (!A.length || !B.length) { statusEl.className = "status err"; statusEl.textContent = `One group is empty (A n=${A.length}, B n=${B.length}). Adjust assignment.`; return false; }
      lastSVGName = `KM_${recs.map(r => r.label).join("_")}_pooled`;
      renderKM_AB(A, B, aStrata.concat(bStrata), `${gLabel}\npooled ${cancers.length} cancers (stratified)${state.months > 0 ? ` · ${state.months}-mo OS` : ""}${subsetTxt ? " · " + subsetTxt : ""}`, xCap);
      lastKM = { A, B };
      const evA = A.reduce((s, x) => s + x.e, 0), evB = B.reduce((s, x) => s + x.e, 0);
      const note = document.createElement("div"); note.className = "sg-legend";
      note.textContent = `${groupLabel("A")}: n=${A.length}, events=${evA} · ${groupLabel("B")}: n=${B.length}, events=${evB}${subsetTxt ? ` · subset: ${subsetTxt}` : ""} · pooled ${cancers.length} cancers (cancer-stratified) · OS · ${splitLabel(state.split)}${state.months > 0 ? ` · capped ${state.months} mo` : ""}`;
      resultEl.appendChild(note);
      return true;
    }

    // ---- 多癌別 per-cancer：每癌別 A vs B 的 HR(log2) 一行 heatmap ----
    function drawHeatmap(cancers, geneVals, recs) {
      const row = []; const pflat = []; let subsetTxt = "";
      cancers.forEach((c, ci) => {
        const { patients, txt } = applySubset(patientsInScope(dataset, [c])); subsetTxt = txt;
        const { A, B } = cohortToAB(buildCohort(patients, geneVals), state.nGenes);
        const evA = A.reduce((s, x) => s + x.e, 0), evB = B.reduce((s, x) => s + x.e, 0);
        let cell;
        if (!A.length || !B.length) cell = { state: "nodata", tip: `${c}: cannot form both groups (A=${A.length}/B=${B.length})` };
        else if (evA + evB === 0) cell = { state: "nodata", tip: `${c}: no events` };
        else {
          const allTm = A.concat(B).map(x => x.tm), allE = A.concat(B).map(x => x.e), allG = A.map(() => 1).concat(B.map(() => 0));
          const lr = logRank(allTm, allE, allG), cox = coxPH1(allTm, allE, allG);
          const log2hr = isFinite(cox.hr) && cox.hr > 0 ? Math.log2(cox.hr) : 0;
          const weak = A.length < 10 || B.length < 10 || (evA + evB) < 10;
          cell = { value: log2hr, state: weak ? "weak" : "ok", stars: "", tip: `${c}: HR=${cox.hr.toFixed(2)} (${cox.ciLow.toFixed(2)}–${cox.ciHigh.toFixed(2)}), n=${A.length}/${B.length}, events=${evA + evB}, p=${lr.p.toPrecision(2)}` };
          if (!weak) pflat.push({ ci, p: lr.p });
        }
        row.push(cell);
      });
      const q = benjaminiHochberg(pflat.map(x => x.p));
      pflat.forEach((x, k) => { const cell = row[x.ci]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      let maxAbs = 0; row.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) maxAbs = Math.max(maxAbs, Math.abs(c.value)); });
      const colorMax = state.colorMax > 0 ? state.colorMax : Math.min(4, Math.max(0.5, Math.ceil(maxAbs * 10) / 10));
      lastSVGName = "HR_heatmap_groups";
      resultEl.innerHTML = heatmapSVG([`${groupLabel("A")} vs ${groupLabel("B")}`], cancers, [row], { colorMax, scheme: state.scheme, legendLabel: "log2 HR", caption: `${recs.map(r => r.label).join("·")}${state.months > 0 ? ` · ${state.months}-mo OS` : ""}${subsetTxt ? "\n" + subsetTxt : ""}` });
      const note = document.createElement("div"); note.className = "sg-legend";
      note.textContent = `${groupLabel("A")} vs ${groupLabel("B")}${subsetTxt ? ` · subset: ${subsetTxt}` : ""} · ★ log-rank FDR q<0.05 ★★<0.01 ★★★<0.001 ★★★★<0.0001 · red=HR>1 (${groupLabel("A")} worse) ${state.scheme === "rb" ? "blue" : "green"}=HR<1 (${groupLabel("A")} better) · grey=no data, faded=too few`;
      resultEl.appendChild(note);
    }

    // ---- 匯出 ----
    // ── Subset screening：多選臨床維度（仿 Clinical Overview）──
    function renderScreenDims() {
      screenDimsBox.innerHTML = "";
      if (!dims.length) { screenDimsBox.innerHTML = '<span class="sg-note">No clinical variables available for this dataset.</span>'; return; }
      dims.forEach(d => {
        const on = state.selectedScreenDims.includes(d.id);
        const lab = document.createElement("label"); lab.className = "sg-dim" + (on ? " on" : "");
        const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = on;
        cb.addEventListener("change", () => {
          if (cb.checked) { if (!state.selectedScreenDims.includes(d.id)) state.selectedScreenDims.push(d.id); }
          else state.selectedScreenDims = state.selectedScreenDims.filter(x => x !== d.id);
          commit(); renderScreenDims();
          scrResultEl.innerHTML = ""; scrStatusEl.className = "status"; scrStatusEl.textContent = "";   // 改選維度 → 清 screening 結果
        });
        lab.appendChild(cb);
        const sp = document.createElement("span"); sp.textContent = d.name; lab.appendChild(sp);
        screenDimsBox.appendChild(lab);
      });
    }

    // ── Subset screening：行=各選取癌別+Pooled，列=各維度兩側，格=該亞組內 A vs B 的 log2 HR ──
    function drawScreening(cancers, geneVals, recs) {
      const cols = [];                                    // 列：選取維度的兩側（兩側都列）
      state.selectedScreenDims.forEach(dimId => {
        const dim = dimById[dimId]; if (!dim) return;
        const asg = assignmentFor(dim); const [bl, al] = labelsOf(dim, asg);
        cols.push({ label: bl, dim, side: "baseline", asg });
        cols.push({ label: al, dim, side: "advanced", asg });
      });
      const rows = cancers.map(c => ({ label: c, cancers: [c], pooled: false }));   // 行：個別癌別
      if (cancers.length > 1) rows.push({ label: "Pooled", cancers: cancers.slice(), pooled: true });   // + Pooled（癌別分層）
      const cells = []; const pflat = [];
      rows.forEach((r, ri) => {
        const rowCells = [];
        cols.forEach((col, ci) => {
          const A = [], B = [], aStr = [], bStr = [];
          r.cancers.forEach(c => {
            const pts = patientsInScope(dataset, [c]).filter(p => classify(col.dim, p.clin[col.dim.field], col.asg) === col.side);
            const { A: a, B: b } = cohortToAB(buildCohort(pts, geneVals), state.nGenes);
            a.forEach(x => { A.push(x); aStr.push(c); });
            b.forEach(x => { B.push(x); bStr.push(c); });
          });
          const evA = A.reduce((s, x) => s + x.e, 0), evB = B.reduce((s, x) => s + x.e, 0);
          let cell;
          if (!A.length || !B.length) cell = { state: "nodata", tip: `${r.label} / ${col.label}: cannot form both groups (A=${A.length}/B=${B.length})` };
          else if (evA + evB === 0) cell = { state: "nodata", tip: `${r.label} / ${col.label}: no events` };
          else {
            const allTm = A.concat(B).map(x => x.tm), allE = A.concat(B).map(x => x.e), allG = A.map(() => 1).concat(B.map(() => 0));
            const useStrat = r.pooled && r.cancers.length > 1;                       // Pooled 行 → 癌別分層
            const cox = useStrat ? coxPH1Stratified(allTm, allE, allG, aStr.concat(bStr)) : coxPH1(allTm, allE, allG);
            const lr = useStrat ? logRankStratified(allTm, allE, allG, aStr.concat(bStr)) : logRank(allTm, allE, allG);
            const log2hr = isFinite(cox.hr) && cox.hr > 0 ? Math.log2(cox.hr) : 0;
            const weak = A.length < 10 || B.length < 10 || (evA + evB) < 10;
            cell = { value: log2hr, state: weak ? "weak" : "ok", stars: "", tip: `${r.label} / ${col.label}: HR=${cox.hr.toFixed(2)}, n=${A.length}/${B.length}, events=${evA + evB}, p=${lr.p.toPrecision(2)}` };
            if (!weak) pflat.push({ ri, ci, p: lr.p });
          }
          rowCells.push(cell);
        });
        cells.push(rowCells);
      });
      const q = benjaminiHochberg(pflat.map(x => x.p));                              // FDR 跨所有有效格子
      pflat.forEach((x, k) => { const cell = cells[x.ri][x.ci]; const st = pStars(q[k]); cell.stars = st === "ns" ? "" : st; cell.tip += `, q=${q[k].toPrecision(2)}`; });
      let maxAbs = 0; cells.forEach(rc => rc.forEach(c => { if (c.state !== "nodata" && isFinite(c.value)) maxAbs = Math.max(maxAbs, Math.abs(c.value)); }));
      const colorMax = state.scrColorMax > 0 ? state.scrColorMax : Math.min(4, Math.max(0.5, Math.ceil(maxAbs * 10) / 10));
      const gLabel = recs.map(r => r.label).join("·");
      scrLastSVGName = "screening_" + recs.map(r => r.label).join("_");
      scrResultEl.innerHTML = heatmapSVG(rows.map(r => r.label), cols.map(c => c.label), cells, { colorMax, scheme: state.scrScheme, legendLabel: "log2 HR", caption: `${gLabel}${state.months > 0 ? ` · ${state.months}-mo OS` : ""}\n${groupLabel("A")} vs ${groupLabel("B")} across subgroups` });
      const note = document.createElement("div"); note.className = "sg-legend";
      note.textContent = `${groupLabel("A")} vs ${groupLabel("B")} · each cell = log2 HR within that subgroup · Pooled row = cancer-stratified · ★ log-rank FDR q<0.05 ★★<0.01 ★★★<0.001 ★★★★<0.0001 · red=HR>1 (${groupLabel("A")} worse) ${state.scrScheme === "rb" ? "blue" : "green"}=HR<1 · grey=no data, faded=too few · exploratory, not an interaction test`;
      scrResultEl.appendChild(note);
      $("#sg-scr-svg").style.display = ""; $("#sg-scr-png").style.display = "";
    }

    async function runScreening() {
      scrResultEl.innerHTML = ""; $("#sg-scr-svg").style.display = "none"; $("#sg-scr-png").style.display = "none";
      const raw = state.genes.slice(0, state.nGenes);
      if (raw.length < state.nGenes || raw.some(g => !g)) { scrStatusEl.className = "status err"; scrStatusEl.textContent = `Enter ${state.nGenes} gene(s) above.`; return; }
      const recs = []; const bad = [];
      raw.forEach(g => { const r = dataset.resolveGene(g); if (r.error || r.multiple) bad.push(g); else recs.push({ rec: r.rec, label: r.rec.symbol || r.rec.gene_id }); });
      if (bad.length) { scrStatusEl.className = "status err"; scrStatusEl.textContent = "Unrecognized gene(s): " + bad.join(", "); return; }
      if (!state.selectedCancers.length) { scrStatusEl.className = "status err"; scrStatusEl.textContent = "Select at least one cancer (shared, above)."; return; }
      if (!state.selectedScreenDims.length) { scrStatusEl.className = "status err"; scrStatusEl.textContent = "Pick at least one clinical variable to screen."; return; }
      const aKeys = sigKeys(state.nGenes).filter(k => state.assign[k] === "A"), bKeys = sigKeys(state.nGenes).filter(k => state.assign[k] === "B");
      if (!aKeys.length || !bKeys.length) { scrStatusEl.className = "status err"; scrStatusEl.textContent = "Assign at least one signature to A and one to B (shared, above)."; return; }
      scrStatusEl.className = "status"; scrStatusEl.textContent = `Fetching ${recs.length} gene(s)…`;
      let geneVals;
      try { geneVals = await Promise.all(recs.map(x => dataset.getGeneValues(x.rec))); }
      catch (e) { scrStatusEl.className = "status err"; scrStatusEl.textContent = "Failed to load gene files: " + e.message; return; }
      scrStatusEl.textContent = "Screening…";
      const cancers = state.cancers.filter(c => state.selectedCancers.includes(c));
      drawScreening(cancers, geneVals, recs);
      scrStatusEl.className = "status"; scrStatusEl.textContent = "Done.";
    }

    function exportPrism() {
      if (!lastKM) return;
      const rows = [];
      lastKM.B.forEach(r => rows.push({ m: r.tm, low: r.e, high: "" }));   // B 當第一欄
      lastKM.A.forEach(r => rows.push({ m: r.tm, low: "", high: r.e }));
      rows.sort((a, b) => a.m - b.m);
      const csv = `Month,${groupLabel("B")},${groupLabel("A")}\n` + rows.map(r => `${r.m.toFixed(3)},${r.low},${r.high}`).join("\n") + "\n";
      dl(new Blob([csv], { type: "text/csv" }), `prism_${lastSVGName}.csv`);
    }
    function dl(blob, name) { const u = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = u; a.download = name; a.click(); URL.revokeObjectURL(u); }
    function svgFile(svg, name) { if (!svg) return; dl(new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml;charset=utf-8" }), name + ".svg"); }
    function pngFile(svg, name) {
      if (!svg) return;
      const vb = svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal : null;
      const W = vb && vb.width ? vb.width : svg.clientWidth, H = vb && vb.height ? vb.height : svg.clientHeight;
      const scale = Math.max(2, Math.ceil(2400 / W));                  // 目標寬 ≥2400px → 高解析度輸出
      const xml = new XMLSerializer().serializeToString(svg);
      const img = new Image();
      img.onload = () => { const c = document.createElement("canvas"); c.width = Math.round(W * scale); c.height = Math.round(H * scale); const ctx = c.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, c.width, c.height); ctx.drawImage(img, 0, 0, c.width, c.height); c.toBlob(b => dl(b, name + ".png"), "image/png"); };
      img.onerror = () => alert("PNG export failed (try SVG).");
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(xml);
    }
    function downloadSVG() { svgFile(resultEl.querySelector("svg"), lastSVGName); }
    function downloadPNG() { pngFile(resultEl.querySelector("svg"), lastSVGName); }
    function downloadScrSVG() { svgFile(scrResultEl.querySelector("svg"), scrLastSVGName); }
    function downloadScrPNG() { pngFile(scrResultEl.querySelector("svg"), scrLastSVGName); }
    $("#sg-run").addEventListener("click", run);
    $("#sg-run-screen").addEventListener("click", runScreening);
    $("#sg-scr-all").addEventListener("click", () => { state.selectedScreenDims = dims.map(d => d.id); commit(); renderScreenDims(); scrResultEl.innerHTML = ""; scrStatusEl.className = "status"; scrStatusEl.textContent = ""; });
    $("#sg-scr-clear").addEventListener("click", () => { state.selectedScreenDims = []; commit(); renderScreenDims(); scrResultEl.innerHTML = ""; scrStatusEl.className = "status"; scrStatusEl.textContent = ""; });
    $("#sg-scr-svg").addEventListener("click", downloadScrSVG);
    $("#sg-scr-png").addEventListener("click", downloadScrPNG);
    $("#sg-prism").addEventListener("click", exportPrism);
    $("#sg-svg").addEventListener("click", downloadSVG);
    $("#sg-png").addEventListener("click", downloadPNG);

    renderGenes();
    renderCancers();
    renderGrid();
    updateControls();
    renderScreenDims();
  },
};
