# miRNA 整合交接單（給前端 Claude Code）

> **一句話**：TCGA miRNA 成熟體表現已經以「偽基因」形式上傳到 R2（與 gene 檔同 bucket、同樣本軸、同 value_scale），
> 資料層完成且已驗證。**你的工作只剩前端**：讓成熟體能像基因一樣被 GOI 搜尋與繪圖。
> 完整格式準據見 `R2_DATASET_FORMAT.md`（§4.1 index 形狀、§6 掛接點）。

---

## 1. 現況（已完成，不用你動）

- **已上傳到** `tcga-data` bucket（與基因同一個），`rclone check --one-way` 回報 **0 differences、2,452 matching files**。
- **純新增**：只多了下面這些**新路徑**，沒有覆蓋任何既有 gene 檔 / manifest / gene_index：
  ```
  expr/mir/<MIMAT>.csv.gz    每個成熟體一檔（檔名=MIMAT accession，shard="mir"，約 2,450 個）
  mirna_index.json.gz         成熟體查表（與 gene_index 同形狀）
  mirna_manifest.json         provenance（assay/unit/value_scale/n_mature…）
  ```
- **格式與 gene 檔完全相容**：值檔是逗號單列、`round(RPM * 1000)` 整數、缺值空欄、對齊 `samples.csv` 順序。
  `value_scale` 與 gene 一致（1000），前端共用同一個 `÷ scale`。

## 2. mirna_index 的形狀（為什麼能直接併入 geneIndex）

```json
{ "by_ensembl": {
    "MIMAT0000089": {"gene_id":"MIMAT0000089","symbol":"hsa-miR-31-5p","biotype":"miRNA",
                     "shard":"mir","file":"MIMAT0000089.csv.gz","assay":"mirna_rpm","mimat":"MIMAT0000089"} },
  "by_symbol": { "HSA-MIR-31-5P": ["MIMAT0000089"] } }
```
- `by_ensembl` key 用 MIMAT；`by_symbol` key 用**大寫**成熟體名（配合 `resolveGene` 把輸入轉大寫）。
- 形狀與 `gene_index` 相同，所以併入後 `resolveGene` / `getGeneValues` **原樣可用**，路徑模板 `expr/{shard}/{file}` 對 `shard="mir"` 也成立。

## 3. 你要改的（就兩處核心 + 各模組微調）

### 3.1 `js/core/dataset.js` — `load()` 末端：併入 miRNA 查表
在載入 `this.geneIndex` 那行**之後**加：
```js
// 併入 miRNA 查表：成熟體與基因走同一條 resolveGene / getGeneValues
try {
  const mi = JSON.parse(await fetchGzipText(this.base + "/mirna_index.json.gz"));
  Object.assign(this.geneIndex.by_ensembl, mi.by_ensembl);
  for (const [sym, arr] of Object.entries(mi.by_symbol))
    (this.geneIndex.by_symbol[sym] ||= []).push(...arr);
} catch (e) { /* 此資料集沒有 miRNA（如未來 OSCC 尚未建），略過 */ }
```

### 3.2 `js/core/dataset.js` — `getGeneValues()`：空欄 → NaN
```js
// 改前：const vals = txt.split(",").map(v => Number(v) / scale);
const vals = txt.split(",").map(v => (v === "" ? NaN : Number(v) / scale));
```
> 基因檔沒有空欄，此改對基因零副作用；對 miRNA 的缺值樣本是**必要**的
> （否則 `Number("")` 會變成 `0`，把「沒測 miRNA」誤當成 RPM=0）。

### 3.3 走 `getGeneValues` 的分析模組（label + 丟 NaN）
`groupCompare.js` / `clinicalOverview.js` / `survival.js` / `survivalGroups.js` / `coxRegression.js`，當 GOI 是 miRNA 時：
- **軸/標題標 `RPM`** 而非 TPM：判斷 `rec.assay === "mirna_rpm"`（或 `rec.biotype === "miRNA"`）→ 標 `RPM` / `log2(RPM+1)`。
- **統計前先濾掉 NaN**（少數樣本無 miRNA）：分組取值後丟 `NaN` 再算 mean/median/SD/Mann–Whitney/Cox。
  `immuneCorr.js` 已有 NaN 處理，照其寫法即可。

> 統計本身對 RPM/TPM 一視同仁（log2(x+1)、U-test、HR 不變），**只有顯示單位與 NaN 過濾**要處理。
> 3.1+3.2 改完，miRNA 就在**所有 GOI 模組**同時可被搜到（因為都走同一條取值路徑）；3.3 是收尾。

## 4. 你改完後怎麼驗（本機，不碰線上站）

1. 套用 3.1/3.2，`node --check js/core/dataset.js`。
2. `python -m http.server 8000`，開 localhost（baseUrl 已指向線上 R2）。
3. GOI 框打 `hsa-miR-21-5p` → 應能 resolve → 抓 `expr/mir/hsa-miR-21-5p.csv.gz` → 畫圖。
4. **手算抽查一個值**確認沒有靜默錯位：挑一個樣本，回它的 isoform 檔把該 MIMAT 的 RPM 加總、×1000 四捨五入，對前端顯示值（÷1000 還原）。
5. 確認軸標 `RPM`（3.3 之後）。都 OK 才部署 Pages。

## 5. 注意

- **不要動** gene 檔 / `gene_index.json.gz` / `manifest.json`；miRNA 是獨立新增。
- 之後若重建 miRNA 重傳，用 `rclone copy`（**不是 sync**）；若成熟體集合有增減、要清 `expr/mir/` 的孤兒檔，請 scope 到該資料夾再處理。
- 完整格式標準（含未來 OSCC 新資料集的做法）見 `R2_DATASET_FORMAT.md`。
