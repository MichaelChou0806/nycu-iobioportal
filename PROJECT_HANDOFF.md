# NYCU IOBioPortal — 專案交接文件 (Project Handoff)

> **給接手的 AI / Claude Code 的話**：這是一個已開發到相當成熟的純前端 TCGA 分析工具。讀完這份就能無縫接手。
> **先看 §0 硬規則與 §13 協作方式**——使用者會直接糾錯，期望你直接承認並修正，不要繞圈或過度道歉。
> **動手前務必 §11 的流程**：先看現況（grep/view 確認精確位置）再改，改完一定驗證（`node --check` + smoke test）。

---

## 0. 硬規則（最重要，先看）

| 項目 | 規則 |
|---|---|
| **UI / 圖表 / 匯出 / 錯誤訊息** | **全英文**（投稿需求，不可妥協） |
| **程式碼註解** | **繁體中文** |
| **與使用者溝通** | **繁體中文**（技術術語可用英文） |
| **最終圖乾淨** | 黃旗/警告（人數不足等）只出現在操作介面，**絕不出現在 heatmap/KM 等最終圖上**（可用淡格/斜線表達，但不寫字警告） |
| 使用者自述 | beginner programmer、"vibe coding"，會直接 flag 錯誤 |

---

## 1. 這是什麼

**NYCU IOBioPortal** — 瀏覽器端、免安裝、**純客戶端運算**的 TCGA pan-cancer 基因表現 ＋ 免疫浸潤 ＋ 生存分析工具。打開網頁就能用，**沒有後端伺服器**，所有運算在瀏覽器內完成。

- 線上版：<https://nycu-iobioportal.pages.dev/>（部署在 Cloudflare Pages）
- 使用者：HNSC / OSCC 研究者 @ NYCU
- **最終目標**：把實驗室自己的 **OSCC RNA-seq 資料**接進來，做出可發表的研究分析

**目前有 7 個分析分頁**（app.js `ANALYSES` 順序）：
1. **Clinical Overview** — 多基因 × 多癌別 × 多臨床維度的表現總覽
2. **Survival (KM)** — 單/多基因高低分組的生存（基礎、掃描式）
3. **Advanced Survival** — 多基因組合分組 + 臨床 subset + 亞組 screening 的生存（進階、假設驗證式）
4. **Cox Regression** — univariate / multivariate Cox（基因 + 臨床因子的 HR forest plot；OS/DSS/DFI/PFI endpoint）
5. **Immune Correlation** — 基因表現 × 免疫浸潤相關
6. **Gene Correlation** — GOI 之間兩兩表現相關（scatter / 上三角矩陣 / 跨癌種；含 miRNA）
7. **Group Comparison** — 基因表現在臨床兩組間的差異（最早的原型）

---

## 2. 架構理念（為什麼這樣設計）

**① 重運算離線、前端輕量**
重的資料處理（33 癌別、36,833 基因、上萬樣本）在使用者 PC 離線跑「一次」，切成小檔放 Cloudflare R2（免費靜態儲存）。前端只 fetch 當下需要的那幾個基因再做輕量統計。→ 前端快、零伺服器成本、可無限擴充使用者。

**② Dataset 抽象、與資料集無關**（整個專案的靈魂）
每個分析器操作的是抽象的 `Dataset` 介面（取基因值、臨床、免疫分數），**不綁死 TCGA**。
→ OSCC 只要做成相同格式、在 `datasets.json` 加一條，就能套用**同一套分析邏輯**，不用重寫。

**③ 模組化：core ＋ analyses ＋ 薄殼**
引擎共用；每個分析自成一檔、互不干擾；`app.js` 只負責註冊與切頁（**隱藏而非銷毀**，保留狀態）。

**一句話原則**：在 TCGA（大、公開）上開發並驗證；OSCC（自己的）以「新增一個資料集」的方式接入，**而非另開專案**。

---

## 3. 資料層（已完成、已上傳 R2）

**R2 公開 URL**：`https://pub-960ed97df8134c4a848e08ddb6cc560e.r2.dev`
**Bucket**：`tcga-data`（需設 CORS：AllowedOrigins `["*"]`、Methods GET/HEAD、Headers `["*"]`）

**檔案清單**：
- `manifest.json` — `value_scale=1000`（儲存 `round(TPM*1000)` 整數，**前端要除以 1000**）、33 癌別、36,833 基因、11,370 樣本
- per-gene `.gz`（`expr/{shard}/{gene}.csv.gz`）— 每個基因一檔，值的順序對齊 `samples.csv` 的樣本軸
- `samples.csv.gz` — 欄位 `sample_id, patient_id, cancer, sample_class`（這是對齊用的軸，**永不修改**）
- `gene_index.json.gz` — `by_ensembl` / `by_symbol` 查表
- `clinical_final.csv.gz` — 11,160 病人；有 `OS, OS.time(天), DSS, PFI, node_status, stage` 等；**無 HPV / ENE**
- `immune_scores.csv.gz` — 列=cell_type、值沿樣本軸（未匹配→空白）
- `immune_manifest.json` — cell_types 清單 ＋ method 分組

**miRNA assay（已上傳，偽基因形式）**
- `expr/mir/<MIMAT>.csv.gz` — 每個成熟 miRNA 一檔（檔名=MIMAT accession，對齊 gene 用 ENSG 的慣例），
  shard="mir"；值=`round(RPM*1000)`、缺值空欄、對齊 samples.csv（與 gene 檔位元組相容）
- `mirna_index.json.gz` — 與 gene_index 同形狀（by_ensembl 用 MIMAT、by_symbol 用大寫名）；
  前端 `load()` 併入 geneIndex → 成熟體可像基因一樣被 GOI 搜尋（打 `hsa-miR-21-5p`）
- `mirna_manifest.json` — assay=mirna_rpm、value_scale=1000、file_key=MIMAT
- **前端整合已完成**（2026-06-27）：`dataset.js` `load()` 併入 `mirna_index`、`getGeneValues` 空欄→NaN；`clinicalOverview`/`groupCompare` 統計前丟 NaN + GOI 是 miRNA 時標 RPM；其餘 GOI 模組本來就 `isFinite` 濾值。對線上 R2 抽 `hsa-miR-21-5p`（MIMAT0000076）驗過 resolve/對齊/scale/NaN。通用格式標準 → `R2_DATASET_FORMAT.md`

**免疫資料細節**（TIMER2.0）：119 種 cell type，分 7 法（TIMER 6 / CIBERSORT 22 / CIBERSORT-ABS 22 / QUANTISEQ 11 / MCPCOUNTER 11 / XCELL 39 / EPIC 8）。對齊率 97.5%（11,089/11,370，未匹配多為 normal，TIMER 只有 tumor）。對齊靠 `barcode15 = str(sample_id)[:15]`。

---

## 4. 離線 Pipeline（在使用者本機，已完成）

依序執行（腳本都在專案 outputs，非前端的一部分）：
1. `organize_tcga_pancancer.py` — GDC STAR counts → 每癌種 TPM 矩陣 + 樣本/基因註解
2. `build_clinical_final.R` — TCGA-CDR(存活) + PanCanAtlas followup(分期) → clinical 表
3. `build_r2_dataset.py` — 切 per-gene gz ＋ manifest/samples/gene_index
4. `build_immune.py` — TIMER2.0 免疫分數對齊樣本軸（有 `KEY_MODE` 開關：`"tcga15"` 截斷 / `"exact"` 給 OSCC 用）
5. `rclone copy … r2:tcga-data --s3-no-check-bucket`
   - **`--s3-no-check-bucket` 必加**：R2 token 沒有 ListBuckets 權限，不加會 403

---

## 5. 前端結構

```
js/core/        引擎（共用，與資料集無關）
  data        R2 fetch / 解壓（fetchGzipText 用 DecompressionStream）
  dataset     Dataset 介面：resolveGene, getGeneValues(對齊+/scale+快取), clinical,
              loadImmune()(lazy 載入→Map<cellType, Float64Array>)
  stats       所有統計（見 §7）
  plots       所有 SVG 繪圖（見 §8）
  dimensions  臨床維度定義：loadDimensions, patientsInScope(tumor/去重/排除redaction),
              classify(binary/numericSplit/ordinal逐級別)
  state       狀態與命名清單（saveNamed/loadNamed/listNames，未來接雲端只改這三個）
  gois        跨分析共用基因清單（localStorage + window 'gois-changed' 事件同步）

js/analyses/    各分析（一檔一個，自給自足）
  clinicalOverview  Clinical Overview
  survival          Survival (KM)
  survivalGroups    Advanced Survival
  coxRegression     Cox Regression ← univariate / multivariate Cox（forest plot）← 最近開發
  immuneCorr        Immune Correlation
  correlation       Gene Correlation ← GOI 兩兩相關（矩陣 / scatter / 跨癌種）
  groupCompare      Group Comparison

app.js          薄殼：ANALYSES 陣列註冊、切頁（HIDE 不銷毀）；分頁名稱用 module 的 .name 自動顯示
```

**註冊順序**（app.js ANALYSES）：`[clinicalOverview, survival, survivalGroups, coxRegression, immuneCorr, correlation, groupCompare]`
**config/**：`datasets.json`（資料來源，含 R2 baseUrl）、`dimensions.tcga.json`（12 個臨床維度）

---

## 6. 七個分析模組

### ① Clinical Overview (`clinicalOverview.js`, ~625 行)
多基因 × 多癌別 × 多臨床維度總覽。log2FC / row z-score heatmap、Expanded/Condensed 模式、ordinal 維度的 baseline/advanced 逐級別選擇器、FDR。維度勾選 UI（**Advanced Survival 的 screening 維度選擇就是仿這個**）。**基因清單跨分析共享**。

### ② Survival (`survival.js`, ~358 行)
KM 生存（OS、tumor-only）。基礎、掃描式（單基因 pan-cancer）。
- **三種分組**：median 50/50、tertile 上下⅓、quartile Q1/Q4（後兩者丟中間組，顯示每組 n）
- **追蹤月數 cutoff**（administrative censoring）
- **路由**：單基因＋單癌別→KM 曲線；單基因＋多癌別→Per-cancer HR heatmap 或 Pooled KM（癌別分層）；多基因→HR heatmap
- `updateControls()` 依情境顯示/隱藏 colors/scale/multi-cancer。localStorage：`tcga-tool:survival`

### ③ Advanced Survival (`survivalGroups.js`, ~701 行) ★ 最近重點
進階、假設驗證式（多基因組合分組 + 臨床亞組）。與基礎 Survival 互補（深度 vs 廣度），共用 stats/plots，不重寫。
- **N 基因組合分組**：1–4 基因，`sigKeys(n)` 產生 2ⁿ 個 signature（HH/HL/LH/LL…）；`state.assign` = dict signature→"A"/"B"/"exclude"；UI 有 `# of genes` 下拉 + 動態基因框 + 2ⁿ 個 A/B/— 指派格子；presets（All-High vs All-Low / vs rest）
- **分組標籤 `groupLabel(which)`**（使用者很在意，**不要用 A/B**）：單一 signature→直接顯示（`HH`）；若一組正好是另一組的補集→`others`（如 `HHH vs others`）；否則列出（`HL / LH`）
- **二分法**：median/tertile/quartile（`Split` 下拉），與基礎共用 `splitGroups`
- **臨床 subset（Single comparison 的「Clinical subset」下拉）**：分組前先過濾到某維度一側（keep only，交集）
- **亞組 screening（核心新功能）**：多選臨床維度（仿 Clinical Overview 勾選）→ 2D heatmap，每個亞組**各自獨立**做 A vs B（非巢狀篩選）。回答「基因分組效應是否被臨床狀態修飾」（effect modification），對 OSCC 單一隊列特別關鍵（用亞組 screening 取代多癌別掃描）。**定位為探索性，非 interaction test**（legend 寫明）
  - **2D heatmap 軸**：X=選取維度的**兩側都列**（N stage→N0+N+；再加 T stage→4 欄）；Y=選取的個別癌別 **+ Pooled 行**（Pooled 癌別分層）；每格=該亞組內 A vs B 的 log2 HR
- **多癌別 3 態路由**（Single comparison）：單癌別→KM；多癌別 Per-cancer→1 列 HR heatmap；多癌別 Pooled→癌別分層 KM
- **版面三區**：共享設定（Genes/Cancers/Groups）→ Single comparison（Clinical subset + Run comparison）→ Subset screening（多選維度 + Run screening）。各自獨立 Run、各自結果區
- **markStale**：共享設定（基因/癌別/分組/cutoff）一改，**兩邊結果都清空**並提示 re-run（避免新舊並存）；改 subset 或 screening 維度只清自己那邊
- 匯出 Prism CSV（single）+ SVG/PNG（single 與 screening 都有，PNG 高解析度縮放到 ≥2400px 寬）
- localStorage：`tcga-tool:survivalGroups`

### ④ Immune Correlation (`immuneCorr.js`, ~415 行)
基因表現 × 免疫浸潤相關（Spearman 預設 / Pearson）、tumor-only、**per-cancer 不 pooling**。
- **三件式 cell-type 篩選**：method 快選按鈕 ＋ 關鍵字搜尋（跨方法）＋ 已選 chips
- **通用維度選擇器（gene × cancer × cellType）**：全單→scatter；剛好一維多→1-D（Heatmap/Bar 下拉）；剛好兩維多→heatmap（可 Swap）；三維都多→報錯
- r 的 95% CI 用 Fisher z；FDR 星號；toolbar 有 `Y max`、`Label °`（標籤角度）
- **標題規則**：1-D heatmap 不含 GOI（Y 軸已有），X 軸標真正癌別；1-D bar 不含 method（Y 軸已是 "Spearman r"）；cell 名在 `(METHOD)` 前折行。localStorage：`tcga-tool:immune`

### ⑤ Group Comparison (`groupCompare.js`, ~181 行)
最早的原型。單基因 × 單癌別，比較表現在臨床兩組間的差異（Mann–Whitney U）。scatter + bar、匯出 Prism CSV。
- **目前限制**：只有 2 個寫死的分組（vital: Alive/Dead、node: N−/N+），未用 dimensions.tcga.json
- **與其他模組不重複**（Clinical Overview 看臨床分佈不碰表現；這個看表現在臨床組間差異）。但很原始；**未來可考慮**：擴展成用所有 dimensions（像 Clinical Overview 選維度），或評估是否與 Clinical Overview 的單維度視圖合併。目前能用，使用者實測 OK，暫留

### ⑥ Cox Regression (`coxRegression.js`)
逐因子做 Cox。**癌種 = 多選 chips**（仿 survival/correlation，可拖排序、Sort A–Z/by N、Select all/none）：
- **1 癌種 → forest plot**（`plots.js` 的 `forestSVG`，log-scale HR 軸、HR=1 參考線、三態 ok/weak/nodata、紅 HR>1 綠 HR<1）
- **≥2 癌種 → HR heatmap**（列=factor、欄=癌種、色=log2 HR；uni 與 multi 各一張；FDR **每癌種欄內**；rg/rb 配色選單；Swap rows/cols）。**點任一格 → 右側 drill 面板畫該癌種、當前 model 的完整 forest**（仿 correlation `.cr-views`，heatmap 留左、drill 開右；drill 自帶 SVG/PNG/CSV，檔名 `cox_forest_{model}_{cancer}_{ep}`；heatmap 主存檔鍵不變）。點擊依 `state.swap` 反推座標：不 swap→cancerIdx=data-c；swap→cancerIdx=data-r
- 計算抽成 `uniItems`/`multiItems`（forest 與 drill 共用，輸出與舊單癌種路徑一致）
- **Model 切換（segmented toggle）**：Univariate（每因子各自一個 Cox，用 `coxPH1`，跨因子 FDR）/ Multivariate（所有選到的因子放進**一個** `coxPH`，complete-case，互相校正，adjusted HR；圖例顯示 effective n / events / **EPV 警告**；p = Wald；heatmap 模式每癌種各建一個模型）
- **Endpoint 下拉（資料驅動）**：clinical 有 `<EP>`+`<EP>.time` 才列入（TCGA 有 OS/DSS/DFI/PFI），OSCC 自動沿用
- 基因 = High vs Low（median/tertile/quartile，與其他模組一致）；臨床因子 = advanced vs baseline（label 寫成自我說明 "X vs Y"，不用 A/B）；無資料維度自動灰掉
- **vital 排除**（= 存活事件本身會完全分離）；**recurrence 預設不選**（post-baseline，名稱後 amber `†` 標註 + 最下面說明）
- 共線/separation → `coxPH` 回 error、UI 擋下並提示移除 factor。匯出 CSV/SVG/PNG。localStorage `tcga-tool:cox`
- **目前限制**：用現有 11 個臨床維度（去 vital）；全是二元，尚未做多級類別 dummy。**待擴充**：數百項臨床 factor + 可搜尋 picker（見 `CLINICAL_TABLE_SPEC.md`）
- 曾加 VIF 共線性診斷後**移除**：二元因子 VIF 天生偏溫和、幾乎都顯示 low overlap 反而誤導（教訓：二元因子別用 VIF）

### ⑦ Gene Correlation (`correlation.js`)
GOI 之間兩兩表現相關（Spearman 預設 / Pearson）、tumor-only、**per-cancer 不 pooling**。骨架仿 immuneCorr。
- **路由（GOI 數 × 癌種數）**：2 基因×1 癌種 → **scatter**（回歸 + 95% 信賴帶 + 上/右邊緣 histogram）；2 基因×多癌種 → **1-D**（r across cancers，heatmap/bar）；>2 基因×1 癌種 → **上三角相關矩陣**（`corrMatrixSVG`，下三角+對角留白；**點上三角格子 → 該對 scatter** 顯示在下方）；>2 基因×多癌種 → 提示收一軸（高維無法一圖）
- **pairwise-complete**：每一對各自丟 NaN 取交集（miRNA 缺值樣本）→ 每格 n 可能不同
- **色階固定 ±1，可手動調**（`Color max`；TCGA r 普遍低，調 0.3/0.5 放大對比）；scheme rb（r>0 紅、r<0 藍）
- **含 miRNA**：偽基因自動可入；Spearman rank-based 讓 RPM×TPM 混合也合理；scatter 軸各標 RPM/TPM。FDR（矩陣跨所有對、1-D 跨癌種）。匯出 SVG/PNG。localStorage `tcga-tool:corr`
- **未來**：Advanced Correlation（臨床 subset 後的相關，如 Advanced Survival 那樣）— 使用者規劃中

---

## 7. 統計方法 (`stats.js`, ~279 行，全部用獨立計算交叉驗證過)

`mannWhitney`, `median`, `mean`, `sd`, `pStars`, `log2FC`, `benjaminiHochberg`, `zscoreRow`,
`kaplanMeier`, `logRank`, `logRankStratified`,
`coxPH1`（單變量 Cox，Breslow + Newton-Raphson → `{beta,hr,se,ciLow,ciHigh,p}`）, `coxPH1Stratified`（癌別分層）,
`coxPH`（多共變量 Cox，Newton-Raphson + 資訊矩陣反矩陣求 SE → per-covariate `{hr,se,ci,p}` + cov + LR test；奇異回 `{error}`；**p=1 時 == coxPH1**）, `chiSquareP`（chi-square 上尾 p，任意 df，給 LR test）,
`pearsonr`, `spearmanr`（皆回 `{r,p,n}`，p 用 Fisher z；spearman 用平均秩處理 ties）

**驗證過的不變量（改 stats 必須保持）**：
- **Cox 對稱性**：`HR_a × HR_b = 1.00000`；完全相同的分組 → HR=1；分離方向正確
- **log-rank**：O1/E1/V/chiSq/p 與手算一致
- **分層**：單一 stratum == 不分層
- **Group 編碼方向**：所有生存模組用 `allG = A.map(()=>1).concat(B.map(()=>0))`（**A=1, B=0**），所以 **HR = A vs B**（A 相對 B 的風險）。改動或新增生存功能必須沿用此編碼，否則 HR 方向會反

---

## 8. 繪圖層 (`plots.js`, ~319 行)

函數：`scatterSVG, barSVG, heatmapSVG, multiBarSVG, kmCurveSVG, corrScatterSVG, corrBarSVG, forestSVG, corrMatrixSVG`
- `forestSVG`＝橫向 forest：log-scale HR 軸、HR=1 參考線、三態
- `kmCurveSVG`：curves 帶 `times[]` 時在 x 軸下方畫 number-at-risk 表（純加法，gated）
- `corrScatterSVG`：`opts.band`＝95% 回歸信賴帶、`opts.marginals`＝上/右邊緣 histogram（都 gated，immuneCorr 不傳→行為不變）
- `corrMatrixSVG`：上三角相關矩陣，下三角+對角留白；上三角格子帶 `class="corr-cell" data-i/data-j`（呼叫端委派點擊）；色階固定 ±colorMax
- `heatmapSVG`：`opts.clickable`（gated）→ 非 nodata 格 `<rect>` 加 `class="hm-cell" data-r/data-c style="cursor:pointer"`（Cox heatmap 點格 drill 用）；**falsy 時輸出與舊版 byte-identical**（其他 6 個呼叫端不傳→不受影響）

**關鍵技術點（踩過坑、別退回去）**：

1. **SVG 尺寸（最重要）**：所有主圖用內聯 `style="width:100%;max-width:${W}px;height:auto"`，**不要只用 width/height 屬性**。原因：屬性會被外部 CSS 的 `height:auto` 蓋掉，導致小 viewBox 圖被瀏覽器放大 4–5 倍（標題看起來比按鈕大好幾倍）。max-width 內聯 style 優先級高，鎖住「最大不超過原始寬度」
2. **多 SVG gradient id 衝突**：每張 SVG 用唯一 `uid`（否則 legend colorbar 在網頁上空白，PNG 卻正常）
3. **標題折行**：`wrapText` ＋ 支援 `\n` 強制換行（`String(caption).split("\n").flatMap(seg=>wrapText(seg,mc))`）。**kmCurveSVG 也支援多行**（上邊距 T 隨行數自適應 `Math.max(36, 14+lines*16)`）
4. **標籤角度**：`opts.labelAngle`（未傳時預設「長>14字垂直 90°、短斜 45°」）；下邊距用 `sin θ` 算
5. **heatmap 三態格子**：ok / weak（淡＋點）/ nodata（斜線）。**heatmapSVG 支援任意行×列**（screening 的 2D 與 immune 的 2D 都用它）
6. **corrBarSVG Y 軸**：手動 `opts.yMax` 優先；否則自動涵蓋所有 r 與 CI、留 15% 餘裕、`niceCeilR`；`clipPath` 裁切防溢出

---

## 9. 最近幾輪完成的細節（這是離開前的最新狀態）

**Advanced Survival 模組（survivalGroups.js）整個建起來**：
- ✅ 改名 `Survival · Custom Groups` → **`Advanced Survival`**（app.js 用 `.name` 自動顯示，不用改 app.js）
- ✅ N 基因組合分組（2ⁿ signature 指派格子）+ presets + 三種二分法
- ✅ 臨床 subset（keep only，交集）+ 多癌別 3 態路由
- ✅ **亞組 screening**（多選維度 → 2D heatmap，X=亞組兩側、Y=癌別+Pooled、各格獨立 A vs B、FDR、Pooled 癌別分層、定位探索性）
- ✅ 版面三區（共享 / Single comparison / Subset screening，各自 Run）+ **markStale**（共享設定改→清兩邊）
- ✅ screening 的 Select all/Clear + SVG/PNG 高解析度匯出

**標題規則大改（使用者很在意，別退回）**：
- ✅ **不用 "A/B"**：KM 曲線圖例、heatmap 行標籤、legend、Prism CSV 全改用 signature（`HH`/`LL`/`others`/`HL / LH`）
- ✅ **heatmap 標題去掉 "HR —"**（數值是 log2 HR，色條已標，figure legend 會說明）；分組對比移到行標籤
- ✅ subset 維度名**用 `dim.name`（如 "T stage"）對齊 Clinical Overview**，不是 `dim.label`（那是 null）也不是 `dim.id`（t_stage）
- ✅ **標題月數標註**：設了 Follow-up months 才在標題顯示 `{n}-mo OS`（沒設則不顯示）；endpoint OS 移到 legend（未來做 DSS/PFI endpoint 選擇時再考慮放回標題）

**更早幾輪（plots/immune）**：SVG 大小修復（max-width）、免疫標題重構、長標籤垂直、bar 自動/手動 Y scale + clipPath、標籤角度可調。

---

## 10. localStorage keys

- 跨分析基因清單：`tcga-tool:gois`（gois.js，配合 window `gois-changed` 事件）
- `tcga-tool:survival` / `tcga-tool:survivalGroups` / `tcga-tool:immune`

---

## 11. 開發 / 驗證 / 部署流程（請照這個）

- **本地**：ES modules 不能用 `file://`，必須起伺服器：`python -m http.server 8000` → <http://localhost:8000>，改檔後 **Ctrl+F5** 強制刷新
- **改碼方式**：先 grep/view 確認精確位置（行號常因前面改動位移）→ 改 → `node --check <file>` 驗語法 → inline Node smoke test 驗邏輯（檢查 tick 值、clipPath、角度、尺寸約束、HR 對稱、strata 對齊等）→ 回報驗證結果
- **部署**：push 到 Cloudflare Pages；使用者更新 = git clone + 覆蓋新檔
- **R2 上傳**：`rclone copy … r2:tcga-data --s3-no-check-bucket`

---

## 12. 路線圖（接下來要做的）

### 進階生存分析（Advanced Survival 的延伸）
- [x] ~~多基因組合分組 HH vs LL~~（**已完成**，含 subset 與 screening）
- [x] ~~**單 / 多變量 Cox 回歸**（納入臨床共變量）~~（**已完成 v1**：Cox Regression 分頁，univariate + multivariate forest plot、OS/DSS/DFI/PFI endpoint、EPV 警告；新增 `coxPH` 多共變量引擎）。**待擴充**：數百項臨床 factor + 可搜尋 picker（需先整理上傳臨床表，見 `CLINICAL_TABLE_SPEC.md`）、多級類別 dummy coding、（可選）gene × clinical 交互項以嚴格驗證 effect modification
- [ ] **ROC / 時間依賴 ROC**（基因當 predictor）
- [ ] **Logistic regression**

### 相關分析
- [x] ~~**Gene Correlation 分頁**~~（**已完成 2026-06-27**）：GOI 兩兩表現相關，scatter（信賴帶+邊緣分布）/ 上三角矩陣（可點→scatter）/ 跨癌種；per-cancer 不 pooling、pairwise-complete；含 miRNA。新增 `corrMatrixSVG`、`corrScatterSVG` 加 band+marginals
- [ ] **Advanced Correlation**：臨床 subset 後的相關（如 Advanced Survival 那樣先篩臨床再算 r）— 使用者規劃中，稍後做
- [x] ~~**miRNA 成熟體前端整合**~~（**已完成 2026-06-27**）：成熟體以偽基因形式入 GOI（見 §3）

### OSCC 資料接入（最終目標）
- [ ] OSCC RNA-seq 做成相同格式（per-gene、clinical、自跑免疫 deconvolution）
- [ ] `datasets.json` 加一條 ＋ OSCC 專屬 `dimensions`（補 ENE）
- [ ] 免疫：OSCC 跑自己的 deconvolution，`build_immune.py` 設 `KEY_MODE="exact"`（OSCC sample_id 非 barcode、直接對）
- [ ] **地雷**：batch effect（最大風險，記錄 batch、相同 pipeline、考慮 ComBat）；小樣本 power（~100 tumor+30 normal，quartile 後每組 ~25，Cox/KM 不穩）；臨床欄位不同
- 策略：**不要等「全部」OSCC 資料到齊**（永遠不會）。100+30 已足夠先接進來、及早暴露架構不合處

### 其他
- [ ] Group Comparison 是否擴展（用所有 dimensions）或與 Clinical Overview 整合（見 §6⑤）
- [ ] 生存模組若干小 bug（使用者提過，未指明）
- [ ] （未來）雲端命名清單：改 `state.js` 的 saveNamed/loadNamed 接 Cloudflare KV（建議 share-code 而非密碼）

### miRNA 分析
-  [  ] miRNA 整合進 GOI 系統（資料已上傳，見 MIRNA_INTEGRATION_BRIEF.md）

---

## 13. 與使用者協作的方式（請照這個來）

- **迭代式、糾錯驅動**：使用者會在過程中直接抓錯，期望你**直接承認哪裡錯、然後修對**，不繞圈、不過度道歉
- **先看現況再動手**（grep/view 確認精確位置），避免改錯地方；行號常因前面改動位移，改前重新定位
- 改完**一定驗證**並回報結果
- 回覆用**繁體中文、簡潔**；給「放哪個檔、覆寫還是新增、要不要 Ctrl+F5」的明確指示
- 圖表改動使用者常貼截圖——**對照截圖裡的 UI 按鈕文字判斷圖的比例**（他若要縮小才能貼，代表原圖在瀏覽器裡很大）
- 標題/標籤的措辭使用者很講究（不要 A/B、不要在標題重複顯然的東西如 HR、用 dimension 的正式名）——改文字前想清楚，沿用既有規則

---

## 14. 快速啟動檢查清單（接手後想動手時）

1. 起本機 server 開來看，確認 6 個分頁都正常（線上版 <https://nycu-iobioportal.pages.dev/> 可對照；注意線上版尚無 Cox Regression）
2. 確認要做的是路線圖哪一項（Cox 回歸？ROC？OSCC 接入？）
3. 相關檔在 `js/analyses/` 或 `js/core/stats.js`、`plots.js`
4. 新分析 = 寫一個自足的 `js/analyses/xxx.js`（暴露 `{id, name, mount}`）＋ 在 `app.js` ANALYSES 註冊一行
5. 新統計 = 加到 `stats.js`，配 smoke test（保持 §7 不變量，尤其 group 編碼 A=1、Cox 對稱）
6. 任何分析都透過 `Dataset` 介面取資料、臨床切法走 `dimensions.tcga.json`，**不綁 TCGA、不寫死分組**（為 OSCC）
7. 改圖記得 §8 的 SVG max-width 與 gradient uid 兩個雷
