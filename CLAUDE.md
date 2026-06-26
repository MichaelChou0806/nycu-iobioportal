# NYCU IOBioPortal — Claude Code 常駐說明

> 詳細交接見 **PROJECT_HANDOFF.md**（架構、每個模組、不變量、最近進度）。本檔是常駐速查與硬規則，每次動手前都適用。

## 這是什麼
純前端、免安裝的 **TCGA pan-cancer 基因表現 + 免疫 + 生存分析工具**，給 HNSC/OSCC 研究者用。所有運算在瀏覽器內，資料檔放 Cloudflare R2，零後端。最終目標：接入實驗室自己的 OSCC RNA-seq cohort。
線上版：https://nycu-iobioportal.pages.dev/

7 個分析分頁（`js/analyses/`，於 `app.js` ANALYSES 註冊）：
Clinical Overview · Survival (KM) · Advanced Survival · Cox Regression · Immune Correlation · Gene Correlation · Group Comparison

## 不可違反的硬規則
1. **語言**：程式顯示內容（UI 文字、圖表、匯出檔、錯誤訊息）一律**英文**（投稿用）；程式碼**註解用繁體中文**；與使用者**溝通用繁體中文**。
2. **純前端、零後端**：所有計算在瀏覽器；資料放 R2。不要引入伺服器/框架/build step。
3. **切法即資料**：臨床分組規則寫在 `config/dimensions.tcga.json`，**不寫死在 JS**。OSCC 用同格式另開一份維度檔即可。
4. **與資料集無關**：分析透過 `js/core/dataset.js` 的 `Dataset` 介面取資料，**不綁死 TCGA**（為了 OSCC 接入）。
5. **狀態集中**：使用者設定的存/讀只透過 `js/core/state.js`（未來接雲端只改 saveNamed/loadNamed/listNames）。
6. **分析模組可插拔**：每個分析是 `js/analyses/` 下一個自足檔，暴露 `{id, name, mount(container, ctx)}`，於 `app.js` ANALYSES 註冊一行。分頁名稱用 module 的 `.name` 自動顯示。
7. **最終圖乾淨**：人數不足等黃旗只出現在操作介面，**絕不出現在 heatmap/KM 等最終圖上**（可用淡格/斜線表達，不寫字警告）。
8. **不大重構**：沿用既有風格，多檔修改謹慎，不要重新設計。

## 關鍵不變量（改動時別破壞，踩過的坑）
- **SVG 尺寸**（plots.js）：主圖用內聯 `style="width:100%;max-width:${W}px;height:auto"`，**不要只用 width/height 屬性**（會被外部 CSS 蓋掉、圖被放大 4–5 倍）。
- **多 SVG gradient id**：每張 SVG 用唯一 uid（否則 legend colorbar 在網頁上空白）。
- **生存 group 編碼**：`allG = A.map(()=>1).concat(B.map(()=>0))`（**A=1, B=0** → HR = A vs B）。新增生存功能沿用，否則 HR 方向反。
- **Cox 對稱性**（stats.js）：`HR_a × HR_b = 1`；相同分組→HR=1；分層單一 stratum == 不分層。改 stats 要保持。
- **維度顯示名**：用 `dim.name`（如 "T stage"），**不是 `dim.label`（null）也不是 `dim.id`（t_stage）**。
- **分組標籤不用 "A/B"**：生存的分組顯示用 signature（HH/LL/others），不要寫 Group A/B。
- **標題不重複顯然資訊**：heatmap 標題不寫 "HR"（色條已標、figure legend 會說明）。

## 如何執行
ES modules 不能用 `file://`，必須起本機伺服器：
```
python -m http.server 8000   # 開 http://localhost:8000，改檔後 Ctrl+F5
```
R2 需設 CORS（`AllowedOrigins:["*"], Methods:["GET","HEAD"]`），否則瀏覽器擋下 fetch。

## 改碼 / 驗證流程
1. **先 grep/view 確認精確位置**（行號常因前面改動位移）。
2. 改。
3. `node --check <file>` 驗語法。
4. inline Node smoke test 驗邏輯（tick 值、clipPath、HR 對稱、strata 對齊、尺寸約束等）。
5. 回報驗證結果；給使用者「放哪個檔、覆寫、Ctrl+F5」的明確指示（繁中、簡潔）。

## 檔案地圖
```
index.html                     外殼
css/style.css                  全域樣式（分析模組多半自帶 scoped 樣式）
js/app.js                      薄殼：載 datasets、來源下拉、ANALYSES 註冊、切頁(隱藏非銷毀)
js/core/data.js                fetchJson / fetchGzipText(DecompressionStream) / parseCSV / downloadText
js/core/dataset.js             Dataset：resolveGene、getGeneValues(對齊/除scale/快取)、clinical、loadImmune()
js/core/stats.js               mannWhitney/median/sd/pStars/log2FC/benjaminiHochberg/zscoreRow、
                               kaplanMeier/logRank/logRankStratified/coxPH1/coxPH1Stratified、
                               coxPH(多共變量)/chiSquareP、pearsonr/spearmanr
js/core/plots.js               scatterSVG/barSVG/heatmapSVG(三態+星號)/multiBarSVG/kmCurveSVG(含at-risk表)/
                               corrScatterSVG(可+band/marginals)/corrBarSVG/forestSVG/corrMatrixSVG(上三角,可點)
js/core/dimensions.js          loadDimensions、patientsInScope(tumor/去重/排除redaction)、classify(binary/numericSplit/ordinal)
js/core/state.js               saveLast/loadLast、saveNamed/loadNamed/listNames、export/importState、reconcile(版本寬容)
js/core/gois.js                跨分析共用基因清單（localStorage + window 'gois-changed'）
js/analyses/clinicalOverview.js  多GOI×多癌種×多維度總覽（log2FC/z-score heatmap、Expanded/Condensed、FDR）
js/analyses/survival.js          Survival (KM)：單/多基因高低分組生存
js/analyses/survivalGroups.js    Advanced Survival：組合分組 + 臨床 subset + 亞組 screening
js/analyses/coxRegression.js     Cox Regression：univariate/multivariate Cox（基因+臨床因子 HR forest、OS/DSS/DFI/PFI endpoint）
js/analyses/immuneCorr.js        基因 × 免疫浸潤相關（per-cancer、三件式 cell 篩選、維度選擇器）
js/analyses/correlation.js       Gene Correlation：GOI 兩兩相關（scatter/上三角矩陣/跨癌種、per-cancer、含 miRNA）
js/analyses/groupCompare.js      最早原型：單基因×單癌種、表現在臨床兩組差異（Mann–Whitney）
config/datasets.json           資料來源（含 R2 baseUrl）；OSCC 加一筆即可
config/dimensions.tcga.json    12 個臨床維度（切法/級別指派都在此）
```

## 目前狀態與路線圖
- 資料管線、R2、7 個分析模組：**完成**（Cox Regression / Gene Correlation 為本機新加，部署狀態看 git）。
- **Cox 回歸（uni/multivariate）、Gene Correlation、miRNA 成熟體前端整合 已完成**。接下來（詳見 PROJECT_HANDOFF.md §12）：Advanced Correlation（臨床 subset 後相關）、擴充臨床 factor（picker，見 `CLINICAL_TABLE_SPEC.md`）→ ROC → **OSCC cohort 接入（最終目標）**。
- OSCC 接入注意：batch effect（最大風險，ComBat）、小樣本 power（~100+30）、`build_immune.py` 設 `KEY_MODE="exact"`、補 ENE 維度。
