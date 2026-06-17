// =====================================================================
// core/plots.js — SVG 繪圖原件（純函式，回傳 SVG 字串）
// 現在：scatter（散布+中位數）、bar（mean±SD）
// 未來可在這裡加：heatmap、KM 曲線、box plot 等
//
// groups 格式（兩者共用）：
//   scatterSVG: [{ name, vals:[...], med, color }]
//   barSVG:     [{ name, m, sd, n, color }]
// yMax 由呼叫端算好（讓兩張圖共用同一 y 軸，方便比較）
// =====================================================================

const CW = 560, CH = 320, M = { t: 34, r: 20, b: 54, l: 60 };

function yScale(v, yMax) { const ph = CH - M.t - M.b; return M.t + ph - (v / yMax) * ph; }
function groupX(idx, n) { const pw = CW - M.l - M.r; return M.l + pw * (idx + 0.5) / n; }

function axis(yMax) {
  let s = "";
  for (let i = 0; i <= 5; i++) {
    const val = yMax * i / 5, y = yScale(val, yMax);
    s += `<line x1="${M.l}" y1="${y}" x2="${CW - M.r}" y2="${y}" stroke="#eef1f4"/>`;
    s += `<text x="${M.l - 8}" y="${y + 4}" text-anchor="end" font-size="11" fill="#6b7480">${val >= 100 ? val.toFixed(0) : val.toFixed(1)}</text>`;
  }
  s += `<line x1="${M.l}" y1="${M.t}" x2="${M.l}" y2="${CH - M.b}" stroke="#c9ced6"/>`;
  s += `<line x1="${M.l}" y1="${CH - M.b}" x2="${CW - M.r}" y2="${CH - M.b}" stroke="#c9ced6"/>`;
  return s;
}
function frame(inner, caption) {
  return `<svg viewBox="0 0 ${CW} ${CH}" xmlns="http://www.w3.org/2000/svg">
    <text x="${CW / 2}" y="18" text-anchor="middle" font-size="12" fill="#6b7480">${caption}</text>
    ${inner}</svg>`;
}
function xLabel(cx, name, n) {
  return `<text x="${cx}" y="${CH - M.b + 20}" text-anchor="middle" font-size="13" fill="#1f2733" font-weight="600">${name}</text>
    <text x="${cx}" y="${CH - M.b + 36}" text-anchor="middle" font-size="11" fill="#6b7480">n=${n}</text>`;
}

export function scatterSVG(groups, yMax, caption) {
  let g = axis(yMax);
  groups.forEach((grp, gi) => {
    const cx = groupX(gi, groups.length), jw = 42;
    grp.vals.forEach(v => {
      const x = cx + (Math.random() - 0.5) * jw;
      g += `<circle cx="${x.toFixed(1)}" cy="${yScale(v, yMax).toFixed(1)}" r="2.4" fill="${grp.color}" fill-opacity="0.5"/>`;
    });
    const my = yScale(grp.med, yMax);
    g += `<line x1="${cx - 34}" y1="${my}" x2="${cx + 34}" y2="${my}" stroke="#1f2733" stroke-width="2.5"/>`;
    g += xLabel(cx, grp.name, grp.vals.length);
  });
  return frame(g, caption);
}

export function barSVG(groups, yMax, caption) {
  let g = axis(yMax); const bw = 70;
  groups.forEach((grp, gi) => {
    const cx = groupX(gi, groups.length);
    const top = yScale(grp.m, yMax), base = yScale(0, yMax);
    g += `<rect x="${cx - bw / 2}" y="${top}" width="${bw}" height="${base - top}" fill="${grp.color}" fill-opacity="0.75"/>`;
    const eTop = yScale(grp.m + grp.sd, yMax), eBot = yScale(Math.max(0, grp.m - grp.sd), yMax);
    g += `<line x1="${cx}" y1="${eTop}" x2="${cx}" y2="${eBot}" stroke="#1f2733" stroke-width="1.5"/>`;
    g += `<line x1="${cx - 10}" y1="${eTop}" x2="${cx + 10}" y2="${eTop}" stroke="#1f2733" stroke-width="1.5"/>`;
    g += `<line x1="${cx - 10}" y1="${eBot}" x2="${cx + 10}" y2="${eBot}" stroke="#1f2733" stroke-width="1.5"/>`;
    g += xLabel(cx, grp.name, grp.n);
  });
  return frame(g, caption);
}

// =====================================================================
// 交付 B 追加：heatmap（缺資料三態 + 星號）、單基因多組 bar
// =====================================================================

// 紅藍發散配色：v 在 [-vmax, vmax]，0=白、正=紅、負=藍
function divColor(v, vmax, scheme) {
  let t = Math.max(-1, Math.min(1, v / (vmax || 1)));
  if (t >= 0) { const r = 255 + (239 - 255) * t, g = 255 + (68 - 255) * t, b = 255 + (68 - 255) * t; return `rgb(${r | 0},${g | 0},${b | 0})`; } // 正 → 紅
  const u = -t;
  if (scheme === "rg") { const r = 255 + (34 - 255) * u, g = 255 + (197 - 255) * u, b = 255 + (94 - 255) * u; return `rgb(${r | 0},${g | 0},${b | 0})`; } // 負 → 綠（HR favorable）
  const r = 255 + (59 - 255) * u, g = 255 + (130 - 255) * u, b = 255 + (246 - 255) * u; return `rgb(${r | 0},${g | 0},${b | 0})`; // 負 → 藍（預設 log2FC）
}

// heatmap
// rows: [label...]（基因）  cols: [label...]（癌種或分組）
// cells[r][c] = { value:Number|null, state:"ok"|"weak"|"nodata", stars:"", tip:"" }
// opts = { colorMax, legendLabel, caption }
export function heatmapSVG(rows, cols, cells, opts) {
  const uid = "u" + Math.random().toString(36).slice(2, 8);   // 唯一 id，避免多張 SVG 的 defs 衝突
  const cw = 36, ch = 26, R = 84;
  const maxRowLen = Math.max(1, ...rows.map(r => String(r).length));
  const maxColLen = Math.max(1, ...cols.map(c => String(c).length));
  const L = Math.max(100, Math.ceil(maxRowLen * 6.5) + 14);   // 左邊距隨最長列標籤
  const colAngle = opts.labelAngle != null ? opts.labelAngle : (maxColLen > 14 ? 90 : 45);  // 欄標籤角度（預設長垂直短斜），可由 opts 覆寫
  const B = Math.max(60, Math.ceil(maxColLen * 6.0 * Math.sin(colAngle * Math.PI / 180)) + 26);  // 下邊距隨角度的垂直投影
  const W = L + cols.length * cw + R;                                    // 寬度只看內容，不被標題撐寬
  const _mc = Math.max(8, Math.floor((W - 16) / 7.2));
  const titleLines = opts.caption ? String(opts.caption).split("\n").flatMap(seg => wrapText(seg, _mc)) : [];
  const T = Math.max(24, titleLines.length * 15 + 12);                   // 上邊距 = 標題折行後的高度
  const H = T + rows.length * ch + B;
  const vmax = opts.colorMax || 1;
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Segoe UI,Roboto,sans-serif">`;
  s += `<defs><pattern id="hatch-${uid}" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
        <rect width="6" height="6" fill="#eceef1"/><line x1="0" y1="0" x2="0" y2="6" stroke="#cbd0d6" stroke-width="1.5"/></pattern>
        <linearGradient id="cbar-${uid}" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stop-color="${divColor(vmax, vmax, opts.scheme)}"/><stop offset="50%" stop-color="#ffffff"/><stop offset="100%" stop-color="${divColor(-vmax, vmax, opts.scheme)}"/></linearGradient></defs>`;
  titleLines.forEach((ln, i) => s += `<text x="${W / 2}" y="${15 + i * 15}" text-anchor="middle" font-size="13" font-weight="600" fill="#1f2733">${ln}</text>`);

  // 格子
  for (let r = 0; r < rows.length; r++) {
    for (let c = 0; c < cols.length; c++) {
      const cell = cells[r][c] || { state: "nodata" };
      const x = L + c * cw, y = T + r * ch;
      if (cell.state === "nodata") {
        s += `<rect x="${x}" y="${y}" width="${cw}" height="${ch}" fill="url(#hatch-${uid})" stroke="#fff"><title>${cell.tip || "no data"}</title></rect>`;
      } else {
        const op = cell.state === "weak" ? 0.45 : 1;
        s += `<rect x="${x}" y="${y}" width="${cw}" height="${ch}" fill="${divColor(cell.value, vmax, opts.scheme)}" fill-opacity="${op}" stroke="#fff"><title>${cell.tip || ""}</title></rect>`;
        if (cell.state === "weak") s += `<circle cx="${x + cw / 2}" cy="${y + ch / 2}" r="2" fill="#9aa3ad"/>`;
        else if (cell.stars) s += `<text x="${x + cw / 2}" y="${y + ch / 2 + 4}" text-anchor="middle" font-size="12" font-weight="700" fill="#1f2733">${cell.stars}</text>`;
      }
    }
  }
  // 列標籤（基因）
  for (let r = 0; r < rows.length; r++)
    s += `<text x="${L - 6}" y="${T + r * ch + ch / 2 + 4}" text-anchor="end" font-size="12" fill="#1f2733">${rows[r]}</text>`;
  // 欄標籤（短的斜 45°、長的垂直 90°）
  for (let c = 0; c < cols.length; c++) {
    const x = L + c * cw + cw / 2, y = T + rows.length * ch + 8;
    s += `<text x="${x}" y="${y}" transform="rotate(${-colAngle} ${x} ${y})" text-anchor="end" font-size="11" fill="#1f2733">${cols[c]}</text>`;
  }
  // 色階
  const lx = W - R + 22, ly = T, lh = Math.min(rows.length * ch, 140);
  s += `<rect x="${lx}" y="${ly}" width="12" height="${lh}" fill="url(#cbar-${uid})" stroke="#c9ced6"/>`;
  s += `<text x="${lx + 16}" y="${ly + 8}" font-size="10" fill="#6b7480">+${vmax.toFixed(1)}</text>`;
  s += `<text x="${lx + 16}" y="${ly + lh / 2 + 3}" font-size="10" fill="#6b7480">0</text>`;
  s += `<text x="${lx + 16}" y="${ly + lh}" font-size="10" fill="#6b7480">-${vmax.toFixed(1)}</text>`;
  if (opts.legendLabel) s += `<text x="${lx + 6}" y="${ly + lh + 16}" font-size="10" fill="#6b7480">${opts.legendLabel}</text>`;
  s += `</svg>`;
  return s;
}

// 單基因多組 bar（mean±SD；缺資料三態；advanced 帶星號）
// groups[c] = { name, m, sd, state:"ok"|"weak"|"nodata", stars }
// opts = { caption, ylabel }
export function multiBarSVG(groups, opts) {
  const n = groups.length;
  const bw = Math.max(12, Math.min(46, Math.round(560 / Math.max(1, n))));
  const gap = Math.round(bw * 0.5);
  const L = 56, T = 30, B = 96, R = 16;
  const plotW = n * (bw + gap) + gap, plotH = 240;
  const W = L + plotW + R, H = T + plotH + B;
  let yMax = 0;
  groups.forEach(g => { if (g.state !== "nodata") yMax = Math.max(yMax, g.m + (g.err || 0)); });
  if (yMax <= 0) yMax = 1; yMax *= 1.08;
  const yS = v => T + plotH - (v / yMax) * plotH;

  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Segoe UI,Roboto,sans-serif">`;
  if (opts.caption) s += `<text x="${W / 2}" y="16" text-anchor="middle" font-size="13" font-weight="600" fill="#1f2733">${opts.caption}</text>`;
  for (let i = 0; i <= 5; i++) { const val = yMax * i / 5, y = yS(val); s += `<line x1="${L}" y1="${y}" x2="${L + plotW}" y2="${y}" stroke="#eef1f4"/><text x="${L - 6}" y="${y + 4}" text-anchor="end" font-size="10" fill="#6b7480">${val >= 100 ? val.toFixed(0) : val.toFixed(1)}</text>`; }
  s += `<line x1="${L}" y1="${T}" x2="${L}" y2="${T + plotH}" stroke="#c9ced6"/><line x1="${L}" y1="${T + plotH}" x2="${L + plotW}" y2="${T + plotH}" stroke="#c9ced6"/>`;
  if (opts.ylabel) s += `<text x="14" y="${T + plotH / 2}" transform="rotate(-90 14 ${T + plotH / 2})" text-anchor="middle" font-size="11" fill="#6b7480">${opts.ylabel}</text>`;

  groups.forEach((g, i) => {
    const cx = L + gap + i * (bw + gap) + bw / 2;
    const lx = cx, ly = T + plotH + 8;
    s += `<text x="${lx}" y="${ly}" transform="rotate(-45 ${lx} ${ly})" text-anchor="end" font-size="10" fill="#1f2733">${g.name}</text>`;
    if (g.state === "nodata") { return; }   // 無資料：只留標籤、不畫 bar
    const op = g.state === "weak" ? 0.45 : 0.85;
    const top = yS(g.m), base = yS(0);
    s += `<rect x="${cx - bw / 2}" y="${top}" width="${bw}" height="${base - top}" fill="${g.color || '#64748b'}" fill-opacity="${op}"/>`;
    const eT = yS(g.m + (g.err || 0)), eB = yS(Math.max(0, g.m - (g.err || 0)));
    s += `<line x1="${cx}" y1="${eT}" x2="${cx}" y2="${eB}" stroke="#1f2733" stroke-width="1.2"/><line x1="${cx - 5}" y1="${eT}" x2="${cx + 5}" y2="${eT}" stroke="#1f2733" stroke-width="1.2"/>`;
    if (g.state === "ok" && g.stars && g.stars !== "ns") s += `<text x="${cx}" y="${eT - 4}" text-anchor="middle" font-size="11" font-weight="700" fill="#1f2733">${g.stars}</text>`;
  });
  s += `</svg>`;
  return s;
}

// =====================================================================
// Kaplan–Meier 曲線（兩組）。curves=[{label,color,km}]，km 為 kaplanMeier() 輸出。
// opts={pText, hrText, xLabel, caption, xMax}
// =====================================================================
export function kmCurveSVG(curves, opts) {
  const L = 58, T = 36, R = 16, B = 80, plotW = 480, plotH = 300;
  const W = L + plotW + R, H = T + plotH + B;
  const xMax = opts.xMax || Math.max(1, ...curves.map(c => c.km.length ? c.km[c.km.length - 1].t : 1));
  const xS = t => L + (t / xMax) * plotW;
  const yS = sv => T + plotH - sv * plotH;
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Segoe UI,Roboto,sans-serif">`;
  if (opts.caption) s += `<text x="${L + plotW / 2}" y="18" text-anchor="middle" font-size="13" font-weight="600" fill="#1f2733">${opts.caption}</text>`;
  for (let k = 0; k <= 4; k++) { const sv = k / 4, y = yS(sv); s += `<line x1="${L}" y1="${y}" x2="${L + plotW}" y2="${y}" stroke="#eef1f4"/><text x="${L - 6}" y="${y + 4}" text-anchor="end" font-size="10" fill="#6b7480">${sv.toFixed(2)}</text>`; }
  for (let k = 0; k <= 5; k++) { const t = xMax * k / 5, x = xS(t); s += `<line x1="${x}" y1="${T + plotH}" x2="${x}" y2="${T + plotH + 4}" stroke="#c9ced6"/><text x="${x}" y="${T + plotH + 16}" text-anchor="middle" font-size="10" fill="#6b7480">${Math.round(t)}</text>`; }
  s += `<line x1="${L}" y1="${T}" x2="${L}" y2="${T + plotH}" stroke="#c9ced6"/><line x1="${L}" y1="${T + plotH}" x2="${L + plotW}" y2="${T + plotH}" stroke="#c9ced6"/>`;
  s += `<text x="14" y="${T + plotH / 2}" transform="rotate(-90 14 ${T + plotH / 2})" text-anchor="middle" font-size="11" fill="#6b7480">Survival probability</text>`;
  s += `<text x="${L + plotW / 2}" y="${T + plotH + 34}" text-anchor="middle" font-size="11" fill="#6b7480">${opts.xLabel || "Time"}</text>`;
  curves.forEach(c => {
    const km = c.km; if (!km.length) return;
    let d = `M ${xS(km[0].t)} ${yS(km[0].surv)}`;
    for (let i = 1; i < km.length; i++) d += ` H ${xS(km[i].t)} V ${yS(km[i].surv)}`;
    d += ` H ${xS(xMax)}`;
    s += `<path d="${d}" fill="none" stroke="${c.color}" stroke-width="2"/>`;
    km.forEach(p => { if (p.nCensor > 0 && p.t > 0) { const x = xS(p.t), y = yS(p.surv); s += `<line x1="${x}" y1="${y - 4}" x2="${x}" y2="${y + 4}" stroke="${c.color}" stroke-width="1.5"/>`; } });
  });
  let ly = T + 8;
  curves.forEach(c => { s += `<rect x="${L + plotW - 158}" y="${ly - 9}" width="12" height="12" fill="${c.color}"/><text x="${L + plotW - 142}" y="${ly + 1}" font-size="11" fill="#1f2733">${c.label} (n=${c.n})</text>`; ly += 18; });
  const sy = T + plotH - 8;
  if (opts.hrText) s += `<text x="${L + 8}" y="${sy - 16}" font-size="11" fill="#1f2733">${opts.hrText}</text>`;
  if (opts.pText) s += `<text x="${L + 8}" y="${sy}" font-size="11" fill="#1f2733">${opts.pText}</text>`;
  s += `</svg>`;
  return s;
}

// =====================================================================
// 相關散點圖：x=免疫分數, y=基因表現，含最小平方回歸線 + r/p/n 標註
// points=[{x,y}]，opts={xLabel,yLabel,r,pText,n,caption}
// =====================================================================
function fmtAxis(v) { const a = Math.abs(v); return a >= 100 ? v.toFixed(0) : a >= 1 ? v.toFixed(1) : a >= 0.01 ? v.toFixed(3) : v.toExponential(1); }
// 標題折行：超過 maxChars 就斷到下一行（置中、字級不變、不撐寬畫布）
function wrapText(text, maxChars) {
  text = String(text);
  if (text.length <= maxChars) return [text];
  const words = text.split(" "); const lines = []; let cur = "";
  for (const w of words) { if (cur && (cur + " " + w).length > maxChars) { lines.push(cur); cur = w; } else cur = cur ? cur + " " + w : w; }
  if (cur) lines.push(cur); return lines;
}

// 自動取漂亮的 r 軸上界（相關係數最大為 1）
function niceCeilR(x) { const steps = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1]; for (const st of steps) if (x <= st + 1e-9) return st; return 1; }
// 刻度數字格式化（去掉多餘的 0）
function fmtTick(v) { return Math.abs(v) < 1e-9 ? "0" : v.toFixed(2).replace(/0+$/, "").replace(/\.$/, ""); }

export function corrScatterSVG(points, opts) {
  const L = 64, R = 18, B = 64, plotW = 460, plotH = 340;
  const W = L + plotW + R;
  const _mc = Math.max(8, Math.floor((W - 16) / 7.5));
  const titleLines = opts.title ? String(opts.title).split("\n").flatMap(seg => wrapText(seg, _mc)) : [];
  const T = Math.max(40, titleLines.length * 18 + (opts.subtitle ? 18 : 0) + 12);
  const H = T + plotH + B;
  const xs = points.map(p => p.x), ys = points.map(p => p.y);
  const xMin = Math.min(...xs), xMax = Math.max(...xs), yMin = Math.min(...ys), yMax = Math.max(...ys);
  const xPad = (xMax - xMin) * 0.05 || 1, yPad = (yMax - yMin) * 0.05 || 1;
  const x0 = xMin - xPad, x1 = xMax + xPad, y0 = yMin - yPad, y1 = yMax + yPad;
  const sx = v => L + (v - x0) / (x1 - x0) * plotW;
  const sy = v => T + plotH - (v - y0) / (y1 - y0) * plotH;
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Segoe UI,Roboto,sans-serif">`;
  titleLines.forEach((ln, i) => s += `<text x="${W / 2}" y="${20 + i * 18}" text-anchor="middle" font-size="14" font-weight="600" fill="#1f2733">${ln}</text>`);
  if (opts.subtitle) s += `<text x="${W / 2}" y="${20 + titleLines.length * 18}" text-anchor="middle" font-size="12" fill="#6b7480">${opts.subtitle}</text>`;
  s += `<line x1="${L}" y1="${T}" x2="${L}" y2="${T + plotH}" stroke="#c9ced6"/><line x1="${L}" y1="${T + plotH}" x2="${L + plotW}" y2="${T + plotH}" stroke="#c9ced6"/>`;
  for (let k = 0; k <= 5; k++) {
    const yv = y0 + (y1 - y0) * k / 5, y = sy(yv);
    s += `<line x1="${L - 4}" y1="${y}" x2="${L}" y2="${y}" stroke="#c9ced6"/><text x="${L - 7}" y="${y + 4}" text-anchor="end" font-size="10" fill="#6b7480">${fmtAxis(yv)}</text>`;
    const xv = x0 + (x1 - x0) * k / 5, x = sx(xv);
    s += `<line x1="${x}" y1="${T + plotH}" x2="${x}" y2="${T + plotH + 4}" stroke="#c9ced6"/><text x="${x}" y="${T + plotH + 16}" text-anchor="middle" font-size="10" fill="#6b7480">${fmtAxis(xv)}</text>`;
  }
  s += `<text x="18" y="${T + plotH / 2}" transform="rotate(-90 18 ${T + plotH / 2})" text-anchor="middle" font-size="11" fill="#6b7480">${opts.yLabel || "Gene expression"}</text>`;
  s += `<text x="${L + plotW / 2}" y="${T + plotH + 34}" text-anchor="middle" font-size="11" fill="#6b7480">${opts.xLabel || "Immune score"}</text>`;
  points.forEach(p => { s += `<circle cx="${sx(p.x).toFixed(1)}" cy="${sy(p.y).toFixed(1)}" r="2.6" fill="#3b82f6" fill-opacity="0.45"/>`; });
  const n = points.length; let mx = 0, my = 0; for (const p of points) { mx += p.x; my += p.y; } mx /= n; my /= n;
  let num = 0, den = 0; for (const p of points) { num += (p.x - mx) * (p.y - my); den += (p.x - mx) ** 2; }
  if (den > 0) { const slope = num / den, intercept = my - slope * mx;
    s += `<line x1="${sx(x0)}" y1="${sy(slope * x0 + intercept)}" x2="${sx(x1)}" y2="${sy(slope * x1 + intercept)}" stroke="#ef4444" stroke-width="2"/>`;
  }
  s += `<text x="${L + plotW - 6}" y="${T + 14}" text-anchor="end" font-size="11" fill="#1f2733">r = ${opts.r.toFixed(3)},  p ${opts.pText},  n = ${opts.n}</text>`;
  s += `</svg>`;
  return s;
}

// =====================================================================
// r 值 bar（0 為中線，正向上負向下）+ 95% CI 誤差棒 + 星號 + 三態
// items=[{label,r,ciLow,ciHigh,stars,state,tip}]，opts={caption,yLabel,scheme}
// =====================================================================
export function corrBarSVG(items, opts) {
  const L = 56, R = 16, bw = 38, gap = 16;
  const plotW = items.length * (bw + gap) + gap, plotH = 300;
  const maxLabelLen = Math.max(1, ...items.map(it => String(it.label).length));
  const labAngle = opts.labelAngle != null ? opts.labelAngle : (maxLabelLen > 14 ? 90 : 40);  // 標籤角度（預設長垂直短斜），可由 opts 覆寫
  const B = Math.max(70, Math.ceil(maxLabelLen * 6.0 * Math.sin(labAngle * Math.PI / 180)) + 26);  // 下邊距隨角度的垂直投影
  const W = L + plotW + R;
  const _mc = Math.max(8, Math.floor((W - 16) / 7.2));
  const titleLines = opts.caption ? String(opts.caption).split("\n").flatMap(seg => wrapText(seg, _mc)) : [];
  const T = Math.max(30, titleLines.length * 15 + 12), H = T + plotH + B;
  // Y 軸範圍：手動 yMax 優先；否則自動涵蓋所有 r 與 CI、再留 15% 餘裕
  let yMax = opts.yMax > 0 ? opts.yMax : 0;
  if (!yMax) {
    const vs = [];
    items.forEach(it => { if (it.state !== "nodata") { vs.push(Math.abs(it.r)); if (isFinite(it.ciLow)) vs.push(Math.abs(it.ciLow)); if (isFinite(it.ciHigh)) vs.push(Math.abs(it.ciHigh)); } });
    yMax = niceCeilR((vs.length ? Math.max(...vs) : 1) * 1.15);
  }
  const sy = r => T + plotH / 2 - r * (plotH / 2 / yMax);
  const neg = opts.scheme === "rg" ? "#22c55e" : "#3b82f6";
  const uid = "b" + Math.random().toString(36).slice(2, 8);
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto" xmlns="http://www.w3.org/2000/svg" font-family="-apple-system,Segoe UI,Roboto,sans-serif">`;
  s += `<defs><clipPath id="${uid}"><rect x="${L}" y="${T}" width="${plotW}" height="${plotH}"/></clipPath></defs>`;  // 資料只畫在繪圖區內，範圍設太小也不溢出
  titleLines.forEach((ln, i) => s += `<text x="${W / 2}" y="${15 + i * 15}" text-anchor="middle" font-size="13" font-weight="600" fill="#1f2733">${ln}</text>`);
  for (const f of [1, 0.5, 0, -0.5, -1]) { const r = f * yMax, y = sy(r); s += `<line x1="${L}" y1="${y}" x2="${L + plotW}" y2="${y}" stroke="${f === 0 ? '#9aa3af' : '#eef1f4'}"/><text x="${L - 6}" y="${y + 4}" text-anchor="end" font-size="10" fill="#6b7480">${fmtTick(r)}</text>`; }
  s += `<text x="16" y="${T + plotH / 2}" transform="rotate(-90 16 ${T + plotH / 2})" text-anchor="middle" font-size="11" fill="#6b7480">${opts.yLabel || "Correlation (r)"}</text>`;
  s += `<g clip-path="url(#${uid})">`;                                   // 以下 bar / 誤差線 / 星號都裁切在繪圖區內
  items.forEach((it, i) => {
    const cx = L + gap + i * (bw + gap) + bw / 2;
    if (it.state === "nodata") {
      s += `<rect x="${cx - bw / 2}" y="${sy(0) - 1}" width="${bw}" height="2" fill="#cbd5e1"/><text x="${cx}" y="${sy(0) - 6}" text-anchor="middle" font-size="9" fill="#9aa3af">n/a</text>`;
    } else {
      const r = it.r, yTop = sy(Math.max(0, r)), yBot = sy(Math.min(0, r)), op = it.state === "weak" ? 0.4 : 1;
      const color = r >= 0 ? "#ef4444" : neg;
      s += `<rect x="${cx - bw / 2}" y="${yTop}" width="${bw}" height="${Math.max(1, yBot - yTop)}" fill="${color}" fill-opacity="${op}"><title>${it.tip || ''}</title></rect>`;
      if (isFinite(it.ciLow) && isFinite(it.ciHigh)) {
        const yL = sy(it.ciHigh), yH = sy(it.ciLow);
        s += `<line x1="${cx}" y1="${yL}" x2="${cx}" y2="${yH}" stroke="#475569" stroke-width="1.2"/><line x1="${cx - 4}" y1="${yL}" x2="${cx + 4}" y2="${yL}" stroke="#475569" stroke-width="1.2"/><line x1="${cx - 4}" y1="${yH}" x2="${cx + 4}" y2="${yH}" stroke="#475569" stroke-width="1.2"/>`;
      }
      if (it.stars) { const yStar = (r >= 0 ? sy(it.ciHigh) : sy(it.ciLow)) - 4; s += `<text x="${cx}" y="${yStar}" text-anchor="middle" font-size="11" fill="#1f2733">${it.stars}</text>`; }
    }
  });
  s += `</g>`;                                                          // 關閉資料 clip
  items.forEach((it, i) => {                                            // X 標籤畫在繪圖區下方、不裁切；角度可調
    const cx = L + gap + i * (bw + gap) + bw / 2;
    s += `<text x="${cx}" y="${T + plotH + 14}" text-anchor="end" font-size="10" fill="#6b7480" transform="rotate(${-labAngle} ${cx} ${T + plotH + 14})">${it.label}</text>`;
  });
  s += `</svg>`;
  return s;
}
