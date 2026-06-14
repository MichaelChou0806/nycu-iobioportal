// =====================================================================
// core/gois.js — 跨分析共享的 GOIs（基因清單）
// ---------------------------------------------------------------------
// 讓使用者只輸入一次基因，clinicalOverview 與 survival 共用同一份。
// 任一頁改動 → 寫入 localStorage 並廣播 "gois-changed"，另一頁監聽後同步。
// =====================================================================

const KEY = "tcga-tool:gois";

export function getGOIs() {
  try { const v = JSON.parse(localStorage.getItem(KEY)); return Array.isArray(v) ? v : []; }
  catch (e) { return []; }
}

export function setGOIs(arr) {
  const clean = [...new Set((arr || []).map(s => String(s).trim()).filter(Boolean))];
  try { localStorage.setItem(KEY, JSON.stringify(clean)); } catch (e) {}
  window.dispatchEvent(new CustomEvent("gois-changed", { detail: clean }));
  return clean;
}

// 訂閱變更（回傳取消訂閱函式）
export function onGOIsChanged(handler) {
  const fn = e => handler(e.detail || getGOIs());
  window.addEventListener("gois-changed", fn);
  return () => window.removeEventListener("gois-changed", fn);
}

// 從輸入框文字解析成基因陣列
export function parseGenes(text) {
  return [...new Set(String(text).split(/[\s,;]+/).map(s => s.trim()).filter(Boolean))];
}
