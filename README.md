# NYCU IOBioPortal

瀏覽器端、免安裝、純客戶端運算的 TCGA pan-cancer 基因表現 ＋ 免疫 ＋ 生存分析工具。
打開網頁就能用，所有運算在瀏覽器內完成，**沒有後端伺服器**。

線上版：<https://nycu-iobioportal.pages.dev/>

---

## 這是什麼

給實驗室同事用的互動式分析平台。輸入基因（GOI），就能在 33 種 TCGA 癌別上做：

- 多基因 × 多癌別 × 多臨床維度的表現總覽
- Kaplan–Meier 生存分析（單／多基因、單／多癌別、pooled）
- 多基因組合分組的進階生存（HH vs LL ＋ 臨床亞組 screening）
- Univariate ／ Multivariate Cox 回歸（基因 ＋ 臨床因子的 HR forest plot）
- 基因表現 × 免疫浸潤相關性
- 基因 × 基因 的表現相關（兩兩矩陣／scatter，含 miRNA）
- 單基因在臨床兩組間的表現差異

不用裝任何東西、不用寫程式、不用跑 R。

---

## 核心設計理念（為什麼這樣架構）

**1. 重運算離線、前端輕量**
重的資料處理（33 癌別、36,833 基因、上萬樣本）在自己 PC 離線跑「一次」，切成小檔案放 Cloudflare R2（免費靜態儲存）。前端只 fetch 當下需要的那幾個基因，再做輕量統計。→ 前端快、零伺服器成本、可無限擴充使用者。

**2. Dataset 抽象、與資料集無關**
每個分析器操作的是抽象的 `Dataset` 介面（取基因值、臨床、免疫分數），**不綁死 TCGA**。
→ 這是為了接入實驗室自己的 OSCC 資料：OSCC 只要做成相同格式，在 `datasets.json` 加一條就能套用同一套分析，**不用重寫任何分析邏輯**。

**3. 模組化：core ＋ analyses ＋ 薄殼**
引擎（資料／統計／繪圖）共用；每個分析自成一個檔；`app.js` 只負責註冊與切頁。
→ 改一個分析不會動到別的；新增分析只要寫一個檔、註冊一行。

---

## 架構地圖

```
離線 (你的 PC)              R2 (靜態儲存)          前端 (瀏覽器)
──────────────            ──────────           ──────────────
build scripts        →    per-gene .gz     →   fetch 需要的基因
  per-gene TPM            clinical .gz           ↓
  clinical               immune .gz          js/core/    引擎
  immune scores          manifest.json       js/analyses/ 各分析
       │                                      app.js      薄殼
   rclone copy ───────────────┘
```

---

## 前端結構

```
js/core/        引擎（共用）
  data        R2 fetch / 解壓
  dataset     Dataset 介面（取基因值、臨床、免疫）
  stats       所有統計（見下）
  plots       所有 SVG 繪圖（heatmap / bar / scatter / KM …）
  dimensions  臨床維度定義
  state       狀態與命名清單
  gois        跨分析共用的基因清單

js/analyses/    各分析（一檔一個，互不干擾）
  Clinical Overview   多基因 × 多癌別 × 多臨床維度總覽（log2FC / z-score heatmap）
  Survival (KM)       KM 生存（log-rank、Cox HR；單一／多癌別／pooled）
  Advanced Survival   多基因組合分組 ＋ 臨床 subset ＋ 亞組 screening
  Cox Regression      univariate／multivariate Cox（HR forest plot；OS/DSS/DFI/PFI endpoint）
  Immune Correlation  基因表現 × 免疫浸潤相關（Spearman / Pearson，per-cancer）
  Gene Correlation    GOI 兩兩表現相關（scatter／上三角矩陣／跨癌種，含 miRNA）
  Group Comparison    單基因 × 單癌別分組比較（U-test、scatter＋bar）

app.js          薄殼：註冊分析、切換頁籤
```

---

## 統計方法

全部用獨立計算交叉驗證過：
Mann–Whitney U、Kaplan–Meier、log-rank（含 stratified）、Cox PH（單變量／多變量、含 stratified）、
Pearson／Spearman（p 值用 Fisher z 轉換）、Benjamini–Hochberg FDR、chi-square（likelihood-ratio test 用）。

---

## 資料 pipeline（離線，已完成）

1. `organize_tcga_pancancer.py` — 整理 GDC 下載
2. `build_clinical_final.R` — 臨床表（OS / DSS / PFI、stage、node …）
3. `build_r2_dataset.py` — 切 per-gene gz ＋ manifest
4. `build_immune.py` — TIMER2.0 免疫分數對齊樣本軸
5. `rclone copy … r2:tcga-data` — 上傳 R2

---

## 怎麼跑

**本地開發**
```
python -m http.server 8000
# → http://localhost:8000   （改檔後 Ctrl+F5 強制刷新）
```

**部署**：push 到 Cloudflare Pages。

---

## 還沒完成 / 路線圖

### 進階生存分析
- [x] 多基因組合分組（HH vs LL）＋ 臨床亞組 screening（**Advanced Survival** 分頁）
- [x] 單／多變量 Cox 回歸（納入臨床共變量；**Cox Regression** 分頁，HR forest plot、OS/DSS/DFI/PFI endpoint）
- [ ] 數百項臨床 factor ＋ 可搜尋 picker（接入大臨床表後）
- [ ] ROC ／ 時間依賴 ROC（基因當 predictor）
- [ ] Logistic regression

### 相關分析
- [x] 基因 × 基因 相關（**Gene Correlation** 分頁：scatter／上三角矩陣／跨癌種，含 miRNA）
- [ ] Advanced Correlation（臨床 subset 後的相關，如 Advanced Survival 那樣）

### OSCC 資料接入（最終目標）
- [ ] 把實驗室 OSCC RNA-seq 做成相同格式（per-gene、clinical、自跑免疫 deconvolution；臨床表格式見 `CLINICAL_TABLE_SPEC.md`）
- [ ] `datasets.json` 加一條 ＋ OSCC 專屬 dimensions
- [ ] 注意事項：batch effect（最大風險，務必記錄並校正）、小樣本統計 power、臨床欄位不同（補 ENE）

### 其他
- [ ] 生存模組若干小 bug
- [ ] （未來）雲端命名清單：架構已留好，把 `state` 的 saveNamed／loadNamed 換接 Cloudflare KV 即可

---

## 一句話原則

> 在 TCGA（大、公開）上開發並驗證分析；OSCC（自己的）以「新增一個資料集」的方式接入，而非另開專案。
