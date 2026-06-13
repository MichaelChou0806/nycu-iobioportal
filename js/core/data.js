// =====================================================================
// core/data.js — 最低層工具：抓取、解壓、CSV 解析
// 全部是無狀態純函式，不認識任何特定資料來源。
// =====================================================================

// 抓純 JSON（例如 manifest.json）
export async function fetchJson(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error("HTTP " + r.status + " — " + url);
  return r.json();
}

// 抓 .gz 檔並在瀏覽器端解壓成文字（用瀏覽器內建 DecompressionStream，不需套件）
export async function fetchGzipText(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error("HTTP " + r.status + " — " + url);
  if (!("DecompressionStream" in window))
    throw new Error("此瀏覽器不支援 DecompressionStream，請改用較新版 Chrome / Edge / Firefox / Safari");
  const stream = r.body.pipeThrough(new DecompressionStream("gzip"));
  return await new Response(stream).text();
}

// CSV 解析（支援引號內逗號）-> 回傳列陣列（每列是字串陣列）
export function parseCSV(text) {
  const rows = [];
  let row = [], field = "", q = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (q) {
      if (c === '"') { if (text[i + 1] === '"') { field += '"'; i++; } else q = false; }
      else field += c;
    } else {
      if (c === '"') q = true;
      else if (c === ",") { row.push(field); field = ""; }
      else if (c === "\n") { row.push(field); rows.push(row); row = []; field = ""; }
      else if (c === "\r") { /* 略過 */ }
      else field += c;
    }
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  return rows;
}

// CSV -> 物件陣列（用第一列當欄名）
export function parseCSVObjects(text) {
  const rows = parseCSV(text);
  const header = rows[0];
  const out = [];
  for (let i = 1; i < rows.length; i++) {
    if (rows[i].length < 2) continue;
    const o = {};
    header.forEach((h, j) => (o[h] = rows[i][j]));
    out.push(o);
  }
  return { header, rows: out };
}

// 觸發瀏覽器下載一個文字檔（給匯出 CSV 用）
export function downloadText(filename, text, mime = "text/csv") {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}
