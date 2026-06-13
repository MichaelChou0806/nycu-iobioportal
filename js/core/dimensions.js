// =====================================================================
// core/dimensions.js — 臨床維度判讀引擎（共用）
// ---------------------------------------------------------------------
// 讀「維度定義」設定，把病人依某維度分到 baseline / advanced / ignore，
// 並計算 coverage（有資料的癌種/人數）與兩組人數（給勾選介面與計算共用）。
//
// 設計：切法都來自設定（不寫死）。每個 cohort 帶自己的定義檔，
// 未來 OSCC 套同格式即可（相同維度沿用、不同維度各自定義）。
// =====================================================================

import { fetchJson } from "./data.js";

export async function loadDimensions(url) { return fetchJson(url); }

// 取得範圍內的病人：選定癌種、只取 tumor、每病人一個、排除 redacted。
// 回傳 [{patient_id, cancer, clin, idx}]，idx 是該 tumor 樣本在 dataset.samples 的位置（之後算基因值用）。
export function patientsInScope(dataset, cancerCodes) {
  const set = new Set(cancerCodes);
  const seen = new Set();
  const out = [];
  for (let i = 0; i < dataset.samples.length; i++) {
    const s = dataset.samples[i];
    if (!set.has(s.cancer) || s.sample_class !== "tumor") continue;
    if (seen.has(s.patient_id)) continue;
    seen.add(s.patient_id);
    const c = dataset.clinical.get(s.patient_id);
    if (!c) continue;
    if (c.redaction && String(c.redaction).trim() !== "") continue;
    out.push({ patient_id: s.patient_id, cancer: s.cancer, clin: c, idx: i });
  }
  return out;
}

// 有序維度：把原始值（如 T2a、Stage IVA、G3）對到某個級別 id；對不到回 null
function levelOfOrdinal(dim, raw) {
  const v = (raw == null ? "" : String(raw)).trim();
  if (!v) return null;
  const up = v.toUpperCase();
  for (const lv of dim.levels) {
    if (lv.prefix && up.startsWith(lv.prefix.toUpperCase())) return lv.id;
    if (lv.match && lv.match.some(m => up === m.toUpperCase())) return lv.id;
  }
  // 羅馬數字（stage）：去掉 STAGE 後取開頭連續 [IVX]
  const romanLevels = dim.levels.filter(l => l.roman);
  if (romanLevels.length) {
    const r = up.replace("STAGE", "").trim();
    const m = r.match(/^[IVX]+/);
    if (m) { const lv = romanLevels.find(l => l.roman === m[0]); if (lv) return lv.id; }
  }
  return null;
}

function defaultAssign(dim, levelId) {
  const lv = dim.levels.find(l => l.id === levelId);
  return lv ? lv.default : "ignore";
}

// 把一個病人的某維度原始值分到 baseline / advanced / null（null = 不納入或無資料）
// assignment：有序維度的「級別 -> 桶」覆寫（沒給就用定義裡的 default）
export function classify(dim, raw, assignment) {
  if (dim.type === "binary") {
    if (dim.numericSplit) {
      if (raw == null || String(raw).trim() === "") return null;
      const x = Number(raw);
      if (!isFinite(x)) return null;
      const cutoff = (assignment && assignment.cutoff != null) ? assignment.cutoff : dim.numericSplit.cutoff;
      return x < cutoff ? "baseline" : "advanced";
    }
    const v = (raw == null ? "" : String(raw)).trim();
    if (v === "") return null;
    const up = v.toUpperCase();
    if (dim.baseline.values.some(x => x.toUpperCase() === up)) return "baseline";
    if (dim.advanced.values.some(x => x.toUpperCase() === up)) return "advanced";
    return null;
  } else { // ordinal
    const lvId = levelOfOrdinal(dim, raw);
    if (!lvId) return null;
    const bucket = (assignment && assignment[lvId]) || defaultAssign(dim, lvId);
    return bucket === "ignore" ? null : bucket;  // "baseline" | "advanced" | null
  }
}

// 對一個維度，在「範圍內病人」上一次算完：欄位覆蓋、可對映人數、兩組人數、各癌種人數
export function analyzeDimension(dim, patients, assignment) {
  let nField = 0, base = 0, adv = 0;
  const cancersField = new Set();
  const perCancer = {};   // cancer -> {base, adv}
  for (const p of patients) {
    const raw = p.clin[dim.field];
    if (raw != null && String(raw).trim() !== "") { nField++; cancersField.add(p.cancer); }
    const b = classify(dim, raw, assignment);
    if (b) {
      if (b === "baseline") base++; else adv++;
      const pc = (perCancer[p.cancer] || (perCancer[p.cancer] = { base: 0, adv: 0 }));
      pc[b === "baseline" ? "base" : "adv"]++;
    }
  }
  return { nField, nMapped: base + adv, nCancers: cancersField.size, base, adv, perCancer };
}

// 依門檻產生警告訊息（空陣列代表沒問題）
export function warnFor(base, adv, thr) {
  const msgs = [];
  if (base === 0 || adv === 0) { msgs.push("one group empty"); return msgs; }
  if (base < thr.minPerGroup) msgs.push(`baseline only ${base}`);
  if (adv < thr.minPerGroup) msgs.push(`advanced only ${adv}`);
  const hi = Math.max(base, adv), lo = Math.min(base, adv);
  if (hi / lo > thr.maxRatio) msgs.push(`imbalanced ${base}:${adv}`);
  return msgs;
}
