// =====================================================================
// app.js — 薄外殼
// ---------------------------------------------------------------------
// 職責：載入 dataset、資料來源下拉、分析註冊表、掛載分析。
// 重要改動：切換分頁改成「隱藏而非銷毀」，所以你建好的篩選切走再切回都還在。
//          （只有切換『資料來源』時才會重建，因為資料變了。）
// =====================================================================

import { loadDatasetConfigs } from "./core/dataset.js";
import { clinicalOverview } from "./analyses/clinicalOverview.js";
import { survival } from "./analyses/survival.js";
import { survivalGroups } from "./analyses/survivalGroups.js";
import { coxRegression } from "./analyses/coxRegression.js";
import { immuneCorr } from "./analyses/immuneCorr.js";
import { correlation } from "./analyses/correlation.js";
import { groupCompare } from "./analyses/groupCompare.js";

// 分析註冊表：新增分析 = import 進來 + 在這裡加一行
const ANALYSES = [
  clinicalOverview,
  survival,
  survivalGroups,
  coxRegression,
  immuneCorr,
  correlation,
  groupCompare,
];

let datasets = [];
let activeDataset = null;
let currentAnalysis = ANALYSES[0];
const mounted = new Map();   // analysisId -> 該分析的容器 DOM（建一次、留著）

async function main() {
  const statusEl = document.getElementById("appStatus");
  try {
    datasets = await loadDatasetConfigs("config/datasets.json");
  } catch (e) {
    statusEl.className = "status err";
    statusEl.textContent = "Failed to load data (CORS?): " + e.message;
    return;
  }
  activeDataset = datasets[0];

  // 資料來源下拉
  const dsel = document.getElementById("datasetSelect");
  datasets.forEach(d => {
    const o = document.createElement("option");
    o.value = d.id; o.textContent = d.name;
    dsel.appendChild(o);
  });
  dsel.addEventListener("change", async () => {
    const d = datasets.find(x => x.id === dsel.value);
    statusEl.style.display = ""; statusEl.className = "status"; statusEl.textContent = "Loading " + d.name + " …";
    try { await d.load(); } catch (e) {
      statusEl.className = "status err"; statusEl.textContent = "Failed to load: " + e.message; return;
    }
    activeDataset = d;
    showStatus();
    rebuildAll();   // 資料換了，重建所有分析面板
  });

  // 分析選單（分頁）
  const nav = document.getElementById("analysisNav");
  ANALYSES.forEach((a, i) => {
    const b = document.createElement("button");
    b.textContent = a.name; b.className = "tab" + (i === 0 ? " active" : "");
    b.addEventListener("click", () => {
      [...nav.children].forEach(c => c.classList.remove("active"));
      b.classList.add("active");
      currentAnalysis = a;
      showAnalysis(a);
    });
    nav.appendChild(b);
  });

  showStatus();
  showAnalysis(currentAnalysis);
}

function showStatus() {
  // dataset 摘要放到 header 副標（dataset-aware，OSCC 換資料來源自動更新）
  const info = document.getElementById("datasetInfo");
  info.textContent = `${activeDataset.name} · ${activeDataset.samples.length.toLocaleString()} samples · ${activeDataset.clinical.size.toLocaleString()} patients · ${Object.keys(activeDataset.geneIndex.by_ensembl).length.toLocaleString()} genes`;
  const s = document.getElementById("appStatus");
  if (s) s.style.display = "none";   // ready 後隱藏 tab 下方那行（loading/error 時才顯示）
}

// 顯示某分析：把其他的藏起來、顯示這個；第一次顯示才建立（之後留著不銷毀）
function showAnalysis(a) {
  const panel = document.getElementById("analysisPanel");
  mounted.forEach(el => (el.style.display = "none"));
  let el = mounted.get(a.id);
  if (!el) {
    el = document.createElement("div");
    panel.appendChild(el);
    mounted.set(a.id, el);
    a.mount(el, { dataset: activeDataset });   // 建立一次（mount 可能是 async）
  }
  el.style.display = "";
}

// 資料來源變更時才整個重建
function rebuildAll() {
  const panel = document.getElementById("analysisPanel");
  panel.innerHTML = "";
  mounted.clear();
  showAnalysis(currentAnalysis);
}

main();
