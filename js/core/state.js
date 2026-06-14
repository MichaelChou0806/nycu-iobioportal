// =====================================================================
// core/state.js — 狀態存讀層
// ---------------------------------------------------------------------
// 把「使用者的設定狀態」存起來、取回來。重點：
//   - 帶 version，並寬容載入（增減維度/欄位也不會壞）
//   - 現在用瀏覽器 localStorage（零後端）
//   - 未來要接雲端帳號清單，只需替換 saveNamed / loadNamed / listNames
//     這三個函式的內部實作（介面不變），其他程式不用動。
// =====================================================================

const PREFIX = "tcga-tool:";          // localStorage key 前綴
const LAST_KEY = PREFIX + "last";     // 「上次操作」自動記憶

// ---- 自動記憶（上次操作）----
export function saveLast(state) {
  try { localStorage.setItem(LAST_KEY, JSON.stringify(state)); } catch (e) { /* 隱私模式可能失敗，忽略 */ }
}
export function loadLast() {
  try { const s = localStorage.getItem(LAST_KEY); return s ? JSON.parse(s) : null; }
  catch (e) { return null; }
}

// ---- 命名清單（使用者自訂名稱，如 "Epi-import"）----
// 現在存在 localStorage；未來換成雲端只改這三個函式內部。
export function saveNamed(name, state) {
  if (!name) return;
  try { localStorage.setItem(PREFIX + "named:" + name, JSON.stringify(state)); } catch (e) {}
}
export function loadNamed(name) {
  try { const s = localStorage.getItem(PREFIX + "named:" + name); return s ? JSON.parse(s) : null; }
  catch (e) { return null; }
}
export function listNames() {
  const out = [];
  try {
    for (let i = 0; i < localStorage.length; i++) {
      const k = localStorage.key(i);
      if (k && k.startsWith(PREFIX + "named:")) out.push(k.slice((PREFIX + "named:").length));
    }
  } catch (e) {}
  return out.sort();
}
export function deleteNamed(name) {
  try { localStorage.removeItem(PREFIX + "named:" + name); } catch (e) {}
}

// ---- 匯出 / 匯入成檔（跨電腦、可分享）----
export function exportState(state, filename = "tcga-tool-settings.json") {
  const blob = new Blob([JSON.stringify(state, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}
export function importStateFromFile(file) {
  return new Promise((resolve, reject) => {
    const r = new FileReader();
    r.onload = () => { try { resolve(JSON.parse(r.result)); } catch (e) { reject(new Error("Not a valid settings JSON file")); } };
    r.onerror = () => reject(new Error("Failed to read file"));
    r.readAsText(file);
  });
}

// =====================================================================
// 寬容套用：把存下來的狀態安全地套回「目前可用的選項」
// ---------------------------------------------------------------------
// saved   : 載入的舊狀態（可能含已被刪除的維度、或缺少新維度）
// avail   : 目前可用的東西 { cancers:[code...], dimIds:[id...] }
// 規則：保留舊狀態裡仍存在的項目（並維持其順序）、把新出現的項目補在後面、
//       丟掉已不存在的項目。 -> 增減 col 也不會壞。
// =====================================================================
export function reconcile(saved, avail) {
  const out = {
    version: 2,
    genes: (saved && Array.isArray(saved.genes)) ? saved.genes.slice() : [],
    mode: (saved && saved.mode) || "expanded",
    cancers: mergeOrdered(saved && saved.cancers, avail.cancers),
    dimensions: mergeOrdered(saved && saved.dimensions, avail.dimIds),
    selectedCancers: null,   // 由各模組決定預設（見下）
    selectedDims: (saved && Array.isArray(saved.selectedDims)) ? saved.selectedDims.filter(id => avail.dimIds.includes(id)) : avail.dimIds.slice(),
    ordinalAssign: (saved && saved.ordinalAssign && typeof saved.ordinalAssign === "object") ? saved.ordinalAssign : {},
    numericCutoff: (saved && saved.numericCutoff && typeof saved.numericCutoff === "object") ? saved.numericCutoff : {},
  };
  // 勾選的癌種：沿用舊的（過濾掉已不存在的）；沒有就預設全選
  if (saved && Array.isArray(saved.selectedCancers)) {
    out.selectedCancers = saved.selectedCancers.filter(c => avail.cancers.includes(c));
    if (out.selectedCancers.length === 0) out.selectedCancers = avail.cancers.slice();
  } else {
    out.selectedCancers = avail.cancers.slice();
  }
  return out;
}

// 舊順序優先，保留仍存在者，新項目補在最後
function mergeOrdered(savedOrder, availList) {
  const avail = new Set(availList);
  const kept = (Array.isArray(savedOrder) ? savedOrder : []).filter(x => avail.has(x));
  const keptSet = new Set(kept);
  const added = availList.filter(x => !keptSet.has(x));
  return kept.concat(added);
}
