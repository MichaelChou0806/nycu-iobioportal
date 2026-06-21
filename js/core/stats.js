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

// =====================================================================
// 生存分析：Kaplan–Meier、Log-rank、單變量 Cox PH（純前端、確定性）
// 複用本檔既有的 erf()。times[]=時間, events[]=1事件/0censor, group/x[]=0或1
// =====================================================================

function normCDF(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
function chiSqP1(x) { return 2 * (1 - normCDF(Math.sqrt(Math.max(0, x)))); } // df=1 上尾 p

// Kaplan–Meier 階梯點：回傳 [{t, surv, atRisk, nEvent, nCensor}]（含起點 t=0, surv=1）
export function kaplanMeier(times, events) {
  const order = times.map((t, i) => i).sort((a, b) => times[a] - times[b]);
  const n = order.length;
  const pts = [{ t: 0, surv: 1, atRisk: n, nEvent: 0, nCensor: 0 }];
  let surv = 1, atRisk = n, i = 0;
  while (i < n) {
    const t = times[order[i]];
    let d = 0, c = 0;
    while (i < n && times[order[i]] === t) { if (events[order[i]]) d++; else c++; i++; }
    if (d > 0) surv *= (1 - d / atRisk);
    pts.push({ t, surv, atRisk, nEvent: d, nCensor: c });
    atRisk -= (d + c);
  }
  return pts;
}

// Log-rank（Mantel–Cox），兩組 group∈{0,1}。回傳 {chiSq, p, O1, E1, V}
export function logRank(times, events, group) {
  const n = times.length;
  const order = times.map((t, i) => i).sort((a, b) => times[a] - times[b]);
  let O1 = 0, E1 = 0, V = 0;
  let atRisk = n, atRisk1 = group.reduce((s, g) => s + (g === 1 ? 1 : 0), 0);
  let i = 0;
  while (i < n) {
    const t = times[order[i]];
    let d = 0, d1 = 0, leave = 0, leave1 = 0, j = i;
    while (j < n && times[order[j]] === t) {
      const g = group[order[j]] === 1 ? 1 : 0;
      if (events[order[j]]) { d++; if (g) d1++; }
      leave++; if (g) leave1++; j++;
    }
    if (d > 0) {
      E1 += d * atRisk1 / atRisk;
      O1 += d1;
      if (atRisk > 1) V += d * (atRisk1 / atRisk) * (1 - atRisk1 / atRisk) * (atRisk - d) / (atRisk - 1);
    }
    atRisk -= leave; atRisk1 -= leave1; i = j;
  }
  const chiSq = V > 0 ? (O1 - E1) * (O1 - E1) / V : 0;
  return { chiSq, p: V > 0 ? chiSqP1(chiSq) : 1, O1, E1, V };
}

// 單變量 Cox PH（x 為 0/1），Breslow 處理 ties，Newton–Raphson。
// 回傳 {beta, hr, se, ciLow, ciHigh, p}
export function coxPH1(times, events, x) {
  const n = times.length;
  const order = times.map((t, i) => i).sort((a, b) => times[a] - times[b]); // 時間遞增
  function scoreInfo(beta) {
    let U = 0, I = 0, S0 = 0, S1 = 0, S2 = 0, i = n - 1;
    while (i >= 0) {                    // 由最大時間往最小累積 risk set（time >= t）
      const t = times[order[i]];
      let j = i, dSumX = 0, dCount = 0;
      while (j >= 0 && times[order[j]] === t) {
        const xv = x[order[j]], w = Math.exp(beta * xv);
        S0 += w; S1 += w * xv; S2 += w * xv * xv;
        if (events[order[j]]) { dSumX += xv; dCount++; }
        j--;
      }
      if (dCount > 0 && S0 > 0) {       // Breslow：同一時間的事件共用 risk set
        const m = S1 / S0;
        U += dSumX - dCount * m;
        I += dCount * (S2 / S0 - m * m);
      }
      i = j;
    }
    return { U, I };
  }
  let beta = 0;
  for (let it = 0; it < 50; it++) {
    const { U, I } = scoreInfo(beta);
    if (!isFinite(I) || I === 0) break;
    const step = U / I;
    beta += step;
    if (Math.abs(step) < 1e-7) break;
  }
  const { I } = scoreInfo(beta);
  const se = I > 0 ? 1 / Math.sqrt(I) : NaN;
  const hr = Math.exp(beta);
  const z = (se && isFinite(se)) ? beta / se : 0;
  const p = (se && isFinite(se)) ? 2 * (1 - normCDF(Math.abs(z))) : 1;
  return { beta, hr, se, ciLow: Math.exp(beta - 1.96 * se), ciHigh: Math.exp(beta + 1.96 * se), p };
}

// =====================================================================
// 分層（stratified）log-rank 與 Cox —— pan-cancer pooled 用（按癌種分層）
// 在每個 stratum（癌種）內各自建 risk set，再加總，控制癌種基線差異。
// =====================================================================

// 分層 log-rank：strata[] 為每個樣本的分層標籤
export function logRankStratified(times, events, group, strata) {
  const buckets = {};
  for (let i = 0; i < times.length; i++) (buckets[strata[i]] || (buckets[strata[i]] = [])).push(i);
  let O = 0, E = 0, V = 0;
  for (const s in buckets) {
    const idx = buckets[s];
    const r = logRank(idx.map(i => times[i]), idx.map(i => events[i]), idx.map(i => group[i]));
    O += r.O1; E += r.E1; V += r.V;
  }
  const chiSq = V > 0 ? (O - E) * (O - E) / V : 0;
  return { chiSq, p: V > 0 ? chiSqP1(chiSq) : 1, O1: O, E1: E, V };
}

// 分層 Cox（單變量 x∈{0,1}，Breslow），偏似然在每個 stratum 內各自累積 risk set
export function coxPH1Stratified(times, events, x, strata) {
  const buckets = {};
  for (let i = 0; i < times.length; i++) (buckets[strata[i]] || (buckets[strata[i]] = [])).push(i);
  const strataList = Object.values(buckets).map(idx => {
    const t = idx.map(i => times[i]), e = idx.map(i => events[i]), xx = idx.map(i => x[i]);
    const order = t.map((_, i) => i).sort((a, b) => t[a] - t[b]);
    return { t, e, x: xx, order };
  });
  function scoreInfo(beta) {
    let U = 0, I = 0;
    for (const s of strataList) {
      let S0 = 0, S1 = 0, S2 = 0, i = s.order.length - 1;
      while (i >= 0) {
        const tt = s.t[s.order[i]];
        let j = i, dSumX = 0, dCount = 0;
        while (j >= 0 && s.t[s.order[j]] === tt) {
          const xv = s.x[s.order[j]], w = Math.exp(beta * xv);
          S0 += w; S1 += w * xv; S2 += w * xv * xv;
          if (s.e[s.order[j]]) { dSumX += xv; dCount++; }
          j--;
        }
        if (dCount > 0 && S0 > 0) { const m = S1 / S0; U += dSumX - dCount * m; I += dCount * (S2 / S0 - m * m); }
        i = j;
      }
    }
    return { U, I };
  }
  let beta = 0;
  for (let it = 0; it < 50; it++) { const { U, I } = scoreInfo(beta); if (!isFinite(I) || I === 0) break; const step = U / I; beta += step; if (Math.abs(step) < 1e-7) break; }
  const { I } = scoreInfo(beta);
  const se = I > 0 ? 1 / Math.sqrt(I) : NaN;
  const hr = Math.exp(beta), z = (se && isFinite(se)) ? beta / se : 0, p = (se && isFinite(se)) ? 2 * (1 - normCDF(Math.abs(z))) : 1;
  return { beta, hr, se, ciLow: Math.exp(beta - 1.96 * se), ciHigh: Math.exp(beta + 1.96 * se), p };
}

// =====================================================================
// 多共變量 Cox 回歸（Breslow ties + Newton-Raphson）── UniV/MtV 用
// ---------------------------------------------------------------------
// X：n×p 設計矩陣（array of arrays，一列一病人、一行一共變量）。
// 類別變數（3 級以上）請先在外部 dummy/indicator 化成多個 0/1 行再傳進來。
// 回傳每個共變量的 {beta,hr,se,ciLow,ciHigh,z,p} 與整體 LR test。
// 不變量：p=1 時結果與 coxPH1 完全一致；資訊矩陣奇異（共線/separation）回 error。
// =====================================================================

// 對稱矩陣（或一般方陣）反矩陣：Gauss-Jordan + partial pivoting；奇異回 null
function matInverse(Ain) {
  const n = Ain.length;
  const A = Ain.map((r, i) => r.slice().concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
  for (let col = 0; col < n; col++) {
    let piv = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
    if (Math.abs(A[piv][col]) < 1e-12) return null;   // 近奇異
    if (piv !== col) { const t = A[piv]; A[piv] = A[col]; A[col] = t; }
    const d = A[col][col];
    for (let c = 0; c < 2 * n; c++) A[col][c] /= d;
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = A[r][col];
      if (f === 0) continue;
      for (let c = 0; c < 2 * n; c++) A[r][c] -= f * A[col][c];
    }
  }
  return A.map(r => r.slice(n));
}

// log Γ（Lanczos 近似）
function gammln(xx) {
  const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let x = xx, y = xx, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) { y += 1; ser += cof[j] / y; }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}
// 正規化上不完全 Gamma Q(a,x)=1-P(a,x)（series / continued fraction，Numerical Recipes）
function gammaQ(a, x) {
  if (x < 0 || a <= 0) return NaN;
  if (x === 0) return 1;
  if (x < a + 1) {
    let ap = a, sum = 1 / a, del = sum;
    for (let nn = 0; nn < 300; nn++) { ap += 1; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-14) break; }
    return 1 - sum * Math.exp(-x + a * Math.log(x) - gammln(a));
  }
  const FPMIN = 1e-300;
  let b = x + 1 - a, c = 1 / FPMIN, d = 1 / b, h = d;
  for (let i = 1; i <= 300; i++) {
    const an = -i * (i - a);
    b += 2; d = an * d + b; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-14) break;
  }
  return Math.exp(-x + a * Math.log(x) - gammln(a)) * h;
}
// chi-square 上尾 p（自由度 df）；df=1 時與 chiSqP1 一致
export function chiSquareP(x, df) { return gammaQ(df / 2, Math.max(0, x) / 2); }

export function coxPH(times, events, X, opts = {}) {
  const n = times.length;
  const p = (X.length && X[0] && X[0].length) || 0;
  const maxIter = opts.maxIter || 50, tol = opts.tol || 1e-8;
  const nEvents = events.reduce((s, e) => s + (e ? 1 : 0), 0);
  if (p === 0) return { error: "no covariates", converged: false, n, events: nEvents, nCov: 0 };

  const order = times.map((_, i) => i).sort((a, b) => times[a] - times[b]); // 時間遞增

  // 一次掃過：回 score 向量 U、資訊矩陣 Info、log partial likelihood
  // 由最大時間往最小累積 risk set（time >= t），同時間事件用 Breslow 共用 risk set。
  function pass(beta) {
    const U = new Array(p).fill(0);
    const Info = Array.from({ length: p }, () => new Array(p).fill(0));
    let logL = 0, S0 = 0;
    const S1 = new Array(p).fill(0);
    const S2 = Array.from({ length: p }, () => new Array(p).fill(0));
    let i = n - 1;
    while (i >= 0) {
      const t = times[order[i]];
      let j = i, dCount = 0, dSumEta = 0;
      const dSumX = new Array(p).fill(0);
      while (j >= 0 && times[order[j]] === t) {        // 先把此時間點所有列加入 risk set
        const idx = order[j], xr = X[idx];
        let eta = 0; for (let k = 0; k < p; k++) eta += beta[k] * xr[k];
        const w = Math.exp(eta);
        S0 += w;
        for (let k = 0; k < p; k++) { S1[k] += w * xr[k]; const wk = w * xr[k]; const row = S2[k]; for (let l = 0; l < p; l++) row[l] += wk * xr[l]; }
        if (events[idx]) { dCount++; dSumEta += eta; for (let k = 0; k < p; k++) dSumX[k] += xr[k]; }
        j--;
      }
      if (dCount > 0 && S0 > 0) {
        for (let k = 0; k < p; k++) {
          const mk = S1[k] / S0;
          U[k] += dSumX[k] - dCount * mk;
          for (let l = 0; l < p; l++) Info[k][l] += dCount * (S2[k][l] / S0 - (S1[k] / S0) * (S1[l] / S0));
        }
        logL += dSumEta - dCount * Math.log(S0);
      }
      i = j;
    }
    return { U, Info, logL };
  }

  let beta = new Array(p).fill(0), converged = false, iterations = 0, singular = false;
  for (let it = 0; it < maxIter; it++) {
    iterations = it + 1;
    const { U, Info } = pass(beta);
    const inv = matInverse(Info);
    if (!inv) { singular = true; break; }
    let maxStep = 0;
    for (let k = 0; k < p; k++) {
      let s = 0; for (let l = 0; l < p; l++) s += inv[k][l] * U[l];
      beta[k] += s; if (Math.abs(s) > maxStep) maxStep = Math.abs(s);
    }
    if (maxStep < tol) { converged = true; break; }
  }

  const final = pass(beta);
  const cov = matInverse(final.Info);
  if (singular || !cov) return { error: "singular information matrix (check collinearity / separation)", converged: false, n, events: nEvents, nCov: p };

  const loglik0 = pass(new Array(p).fill(0)).logL;
  const hr = [], se = [], ciLow = [], ciHigh = [], z = [], pval = [];
  for (let k = 0; k < p; k++) {
    const v = cov[k][k], s = v > 0 ? Math.sqrt(v) : NaN;
    se[k] = s; hr[k] = Math.exp(beta[k]);
    ciLow[k] = Math.exp(beta[k] - 1.96 * s); ciHigh[k] = Math.exp(beta[k] + 1.96 * s);
    const zz = (s && isFinite(s)) ? beta[k] / s : 0;
    z[k] = zz; pval[k] = (s && isFinite(s)) ? 2 * (1 - normCDF(Math.abs(zz))) : 1;
  }
  const lr = 2 * (final.logL - loglik0);
  return {
    beta: beta.slice(), hr, se, ciLow, ciHigh, z, p: pval, cov,
    loglik: final.logL, loglik0, lr, lrDf: p, lrP: chiSquareP(lr, p),
    n, events: nEvents, nCov: p, iterations, converged,
  };
}

// =====================================================================
// 相關係數：Pearson / Spearman（p 用 Fisher z-transformation，normal 近似）
// 供 GOI 表現 × 免疫分數相關分析。呼叫前請先濾掉任一為 NaN 的配對。
// =====================================================================

// 平均 rank（處理 ties，1-based）
function rankAvg(arr) {
  const idx = arr.map((_, i) => i).sort((a, b) => arr[a] - arr[b]);
  const ranks = new Array(arr.length);
  let i = 0;
  while (i < idx.length) {
    let j = i; while (j < idx.length && arr[idx[j]] === arr[idx[i]]) j++;
    const avg = (i + j - 1) / 2 + 1;
    for (let k = i; k < j; k++) ranks[idx[k]] = avg;
    i = j;
  }
  return ranks;
}
function pearsonCore(x, y) {
  const n = x.length; let sx = 0, sy = 0;
  for (let i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
  const mx = sx / n, my = sy / n; let num = 0, dx = 0, dy = 0;
  for (let i = 0; i < n; i++) { const a = x[i] - mx, b = y[i] - my; num += a * b; dx += a * a; dy += b * b; }
  return (dx > 0 && dy > 0) ? num / Math.sqrt(dx * dy) : 0;
}
// Fisher z：r 的雙尾 p（n<4 或 |r|=1 時退化處理）
function corrP(r, n) {
  if (n < 4) return 1;
  if (Math.abs(r) >= 1) return 0;
  const z = 0.5 * Math.log((1 + r) / (1 - r)), se = 1 / Math.sqrt(n - 3);
  return 2 * (1 - normCDF(Math.abs(z) / se));
}
export function pearsonr(x, y) { const n = x.length; const r = pearsonCore(x, y); return { r, p: corrP(r, n), n }; }
export function spearmanr(x, y) { const n = x.length; const r = pearsonCore(rankAvg(x), rankAvg(y)); return { r, p: corrP(r, n), n }; }
