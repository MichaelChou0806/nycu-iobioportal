// =====================================================================
// core/dataset.js — 資料來源（Dataset）抽象
// ---------------------------------------------------------------------
// 一個 Dataset 物件封裝「某個資料來源」的全部存取方式：
//   載入 manifest / samples / clinical / gene_index、解析基因、抓基因值。
//
// 關鍵設計：分析模組只跟 Dataset「要資料」，不直接寫網址。
// 因為我們的 OSCC cohort 會用同一條 pipeline 產生（格式相同），
// 所以同一個 Dataset 類別就能用在 OSCC，只要在 datasets.json 多加一筆設定。
// =====================================================================

import { fetchJson, fetchGzipText, parseCSVObjects } from "./data.js";

export class Dataset {
  constructor(config) {
    this.config = config;          // { id, name, type, baseUrl }
    this.base = config.baseUrl;
    this.loaded = false;
    this.geneCache = new Map();     // gene_id -> TPM 陣列（避免重複下載）
  }

  get id() { return this.config.id; }
  get name() { return this.config.name; }
  get cancers() { return this.manifest.cancers; }            // [{code,n_tumor,n_normal,has_normal}]
  get clinicalFields() { return this.manifest.clinical_fields; }

  // 載入這個資料來源的基礎檔（只做一次）
  async load() {
    if (this.loaded) return;
    this.manifest = await fetchJson(this.base + "/manifest.json");

    // samples：順序就是每個基因檔裡數值的順序
    const s = parseCSVObjects(await fetchGzipText(this.base + "/" + this.manifest.sample_file));
    this.samples = s.rows;   // [{sample_id, patient_id, cancer, sample_class}]

    // clinical：以 patient_id 建索引
    const c = parseCSVObjects(await fetchGzipText(this.base + "/" + this.manifest.clinical_file));
    this.clinical = new Map();
    c.rows.forEach(r => this.clinical.set(r.patient_id, r));

    // gene_index：symbol / ENSG -> 檔案位置
    this.geneIndex = JSON.parse(await fetchGzipText(this.base + "/" + this.manifest.gene_index_file));

    this.loaded = true;
  }

  // 把使用者輸入（symbol 或 ENSG）解析成基因紀錄
  // 回傳：{rec} 找到 / {multiple:[...]} 重名候選 / {error}
  resolveGene(query) {
    const q = (query || "").trim().toUpperCase();
    if (!q) return { error: "請輸入基因" };
    if (/^ENSG\d+/.test(q)) {
      const base = q.split(".")[0];
      const rec = this.geneIndex.by_ensembl[base];
      return rec ? { rec } : { error: `找不到 ${base}` };
    }
    const bases = this.geneIndex.by_symbol[q];
    if (!bases) return { error: `找不到 symbol「${q}」` };
    if (bases.length === 1) return { rec: this.geneIndex.by_ensembl[bases[0]] };
    return { multiple: bases.map(b => this.geneIndex.by_ensembl[b]) };
  }

  // 抓某基因的 TPM 陣列（順序對齊 this.samples）
  async getGeneValues(rec) {
    if (this.geneCache.has(rec.gene_id)) return this.geneCache.get(rec.gene_id);
    const url = `${this.base}/expr/${rec.shard}/${rec.file}`;
    const txt = await fetchGzipText(url);
    const scale = this.manifest.value_scale || 1;
    const vals = txt.split(",").map(v => Number(v) / scale);
    if (vals.length !== this.samples.length)
      throw new Error(`基因檔長度 ${vals.length} 與樣本數 ${this.samples.length} 不符`);
    this.geneCache.set(rec.gene_id, vals);
    return vals;
  }

  // 載入免疫分數（lazy；只做一次）。對齊 this.samples 軸（與基因值同序）。
  // 回傳 { cellTypes, methods, get(cellType)->Float64Array, has(cellType) }
  async loadImmune() {
    if (this.immune) return this.immune;
    const manifest = await fetchJson(this.base + "/immune_manifest.json");
    const txt = await fetchGzipText(this.base + "/immune_scores.csv.gz");
    const n = this.samples.length;
    const map = new Map();
    txt.split("\n").forEach(line => {
      if (!line) return;
      const parts = line.split(",");           // cell_type, v1, v2, ... vN（空字串=NaN）
      const arr = new Float64Array(n);
      for (let i = 0; i < n; i++) { const v = parts[i + 1]; arr[i] = (v === undefined || v === "") ? NaN : Number(v); }
      map.set(parts[0], arr);
    });
    if (manifest.n_samples !== n) console.warn(`免疫 n_samples ${manifest.n_samples} 與樣本數 ${n} 不符`);
    this.immune = { cellTypes: manifest.cell_types, methods: manifest.methods, get: ct => map.get(ct), has: ct => map.has(ct) };
    return this.immune;
  }
}

// 依 datasets.json 建立所有 Dataset（先載入第一個，其餘等被選到再載入）
export async function loadDatasetConfigs(configUrl) {
  const cfg = await fetchJson(configUrl);
  const list = cfg.datasets.map(d => new Dataset(d));
  await list[0].load();
  return list;
}
