// =====================================================================
// core/stats.js — 統計函式（純函式，無 DOM）
// 現在：Mann–Whitney U + 基本描述統計
// 未來可在這裡加：log-rank、Cox、correlation 等
// =====================================================================

// Mann–Whitney U（含 tie 校正與連續性校正的常態近似，雙尾 p）
export function mannWhitney(a, b) {
  const n1 = a.length, n2 = b.length, N = n1 + n2;
  const comb = a.map(v => ({ v, g: 0 })).concat(b.map(v => ({ v, g: 1 })));
  comb.sort((x, y) => x.v - y.v);
  const ranks = new Array(N);
  let tie = 0, i = 0;
  while (i < N) {
    let j = i;
    while (j + 1 < N && comb[j + 1].v === comb[i].v) j++;
    const avg = (i + j) / 2 + 1;             // 1-based 平均秩
    for (let k = i; k <= j; k++) ranks[k] = avg;
    const t = j - i + 1;
    if (t > 1) tie += t * t * t - t;
    i = j + 1;
  }
  let R1 = 0;
  for (let k = 0; k < N; k++) if (comb[k].g === 0) R1 += ranks[k];
  const U1 = R1 - n1 * (n1 + 1) / 2, U2 = n1 * n2 - U1, U = Math.min(U1, U2);
  const meanU = n1 * n2 / 2;
  const sdU = Math.sqrt((n1 * n2 / 12) * ((N + 1) - tie / (N * (N - 1))));
  let p = 1;
  if (sdU > 0) {
    const zc = (Math.abs(U - meanU) - 0.5) / sdU;   // 連續性校正
    p = 2 * (1 - normCdf(Math.abs(zc)));
    p = Math.min(1, Math.max(0, p));
  }
  return { U, n1, n2, p };
}

function normCdf(x) { return 0.5 * (1 + erf(x / Math.SQRT2)); }
function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}

export function median(arr) {
  const s = [...arr].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
export function mean(arr) { return arr.reduce((a, b) => a + b, 0) / arr.length; }
export function sd(arr) {
  const m = mean(arr);
  return Math.sqrt(arr.reduce((s, v) => s + (v - m) * (v - m), 0) / (arr.length - 1 || 1));
}
export function pStars(p) {
  return p < 1e-4 ? "****" : p < 1e-3 ? "***" : p < 1e-2 ? "**" : p < 0.05 ? "*" : "ns";
}

// =====================================================================
// 交付 B 追加：log2 fold change、Benjamini–Hochberg FDR、row z-score
// =====================================================================

// log2 fold change（advanced 對 baseline，用中位數；加 pseudocount 避免除以 0 / log0）
export function log2FC(baselineVals, advancedVals, pseudo = 0.01) {
  const mb = median(baselineVals) + pseudo;
  const ma = median(advancedVals) + pseudo;
  return Math.log2(ma / mb);
}

// Benjamini–Hochberg：輸入 p 陣列（可含 null/NaN，代表無檢定），回傳同長度的 q（無檢定處回 null）
export function benjaminiHochberg(pvals) {
  const idx = [];
  pvals.forEach((p, i) => { if (p != null && isFinite(p)) idx.push(i); });
  const m = idx.length;
  const q = pvals.map(() => null);
  if (m === 0) return q;
  // 依 p 由小到大排序
  idx.sort((a, b) => pvals[a] - pvals[b]);
  // 由大到小掃，維持單調性
  let prev = 1;
  for (let k = m - 1; k >= 0; k--) {
    const i = idx[k];
    const val = pvals[i] * m / (k + 1);
    prev = Math.min(prev, val);
    q[i] = Math.min(1, prev);
  }
  return q;
}

// 對一列數值做 z-score（忽略 null/NaN；全相同或不足則回全 0/null）
export function zscoreRow(vals) {
  const ok = vals.filter(v => v != null && isFinite(v));
  if (ok.length < 2) return vals.map(v => (v == null || !isFinite(v)) ? null : 0);
  const m = ok.reduce((a, b) => a + b, 0) / ok.length;
  const s = Math.sqrt(ok.reduce((a, b) => a + (b - m) * (b - m), 0) / (ok.length - 1));
  if (s === 0) return vals.map(v => (v == null || !isFinite(v)) ? null : 0);
  return vals.map(v => (v == null || !isFinite(v)) ? null : (v - m) / s);
}
