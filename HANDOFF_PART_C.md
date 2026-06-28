# HANDOFF — Cox Regression 多癌種 heatmap「Part C：點格 → forest drill-down」

> ⚠️ **TRANSIENT（一次性交接檔）。Part C 完成後請刪除本檔**——耐久記錄在 `memory/cox-regression-design.md` 與 `PROJECT_HANDOFF.md §6⑥`。留著會過時/誤導。
> 新 session 開機請先讀：`CLAUDE.md`（硬規則）、`PROJECT_HANDOFF.md`（架構）、memory `cox-regression-design`，再讀本檔。
> 全部只動 `js/analyses/coxRegression.js` 與 `js/core/plots.js` 的 `heatmapSVG`。

---

## 0. 現況（commit hash）
- `9351173` Clinical Overview 加 Normal vs Tumor
- `555c6e4` **Cox Part A**：癌種多選 chips；1 癌種→forest（原狀）、≥2→**Univariate HR heatmap**
- `28a942a` **Cox Part B**：**Multivariate HR heatmap**（每癌種各建一個 coxPH）＋ rg/rb 配色選單＋「Swap」→「Swap rows/cols」
- HEAD = `28a942a`，working tree 乾淨（除本檔）。**尚未 push**（使用者自己 push）。

## 1. Part A/B 已完成（Part C 的基礎，別重做）
`coxRegression.js` 目前：
- **癌種 = 多選 chips**：`state.cancers`（順序）、`state.selectedCancers`、helper `selCancers()`；`renderCancers()`（仿 survival/correlation，可拖、Sort A–Z/by N、Select all/none）。
- **路由（`run()`）**：`selCancers().length===1` → forest（`renderUnivariate`/`renderMultivariate`，已改成吃 `cancer` 參數）；`>=2` → heatmap（model uni→`drawUniHeatmap`、multi→`drawMultiHeatmap`）。
- **`drawUniHeatmap(recs, selDims, geneVals, cancers, epLabel, unknown)`**：`factors = [...genes(kind:"gene"), ...clinical dims(kind:"dim")]`；每 (factor,cancer) 用 `fitRow`→`coxPH1` 算 → `cells[factorIdx][cancerIdx] = {value:log2 HR, state, stars, tip, _hr}`；**FDR 每癌種欄內**；`heatmapSVG`；colorMax auto(cap 4)；`scheme: state.scheme`。
- **`drawMultiHeatmap(...)`**：每癌種 `cancerMultiHR(cancer, factors, geneVals)` → `{byFactor: {factorIdx:{hr,p,weak}}}`（該癌種無資料/單組的 factor 不納入模型→該格 nodata；奇異/cases 不足→該欄全 nodata；EPV<10→weak）；FDR 每欄內；同樣 render。
- **swap**：`state.swap`。轉置：`grid = cancers.map((_,ci)=>factors.map((_,ri)=>cells[ri][ci]))`，並 `rLab=cancers, cLab=facLab`（不 swap 時 `rLab=facLab, cLab=cancers, grid=cells`）。
- **配色**：`state.scheme`（"rg"/"rb"），`#cx-scheme` 下拉（heatmap 時才顯示）。
- **狀態**：`lastView`（"forest"|"heatmap"）、`lastHeatmap = {factors, cancers, cells}`、`lastSVGName`、`lastModel`、`lastCancer`。
- **匯出**：`exportCSV` 在 `lastView==="heatmap"` 時輸出 factor×cancer 的 HR 矩陣（用 `cell._hr`，nodata→NA）；SVG/PNG 走 `firstSVG()`（= `resultEl` 內的 svg）。
- **heatmapSVG 目前不可點**（cell 是純 `<rect>`，無 data 屬性）。

驗證已過：`node --check`；R2 實測 univariate HR heatmap 結構/swap/三態/星號；multivariate per-cancer 路徑（HNSC n=520，KDELR1/2 adj HR 1.39/1.43、converged）。

## 2. Part C 要做的功能
**點 heatmap 的格子 → 開「該癌種的 forest plot」**，且：
- (a) **左右並排**：heatmap 留在左、被點的 forest 開在右側獨立面板（heatmap 不消失，方便連點別格查看）。**仿 `correlation.js` 的 `.cr-views` flex + matrix→scatter 模式**（那邊已驗過好用）。
- (b) **forest 面板有自己的存圖鍵**，且 UI **清楚標示是 heatmap 還是 forest**（例：面板標題「Forest · {cancer}」，存檔檔名含 cancer + "forest"；heatmap 的主存檔鍵維持存 heatmap）。
- (c) forest 內容 = **該癌種、當前 model（uni 或 multi）** 的完整 forest（所有 factor）。
- (d) 容易回到 heatmap（heatmap 一直在；面板可被下一次點擊取代）。

## 3. 預計修改的檔案與函式
### `js/core/plots.js` — `heatmapSVG(rows, cols, cells, opts)`
- 加 **gated** `opts.clickable`（或 `opts.cellAttrs`）：為**非 nodata** 的 cell `<rect>` 加 `class="hm-cell" data-r="${r}" data-c="${c}" style="cursor:pointer"`。
- **不可破壞**：`opts.clickable` falsy 時輸出**與現在 byte-identical**（其他 6 個呼叫端：survival HR heatmap、advanced survival screening、immune 1D/2D、clinicalOverview、correlation 1D、cox uni/multi 都不傳此 opt）。

### `js/analyses/coxRegression.js`
- HTML：把 `#cx-result` 包進一個 flex 容器，並加 `#cx-drill` 面板（仿 correlation `.cr-views > #cr-result + #cr-pair`，空面板 `:empty` 隱藏）。
- `drawUniHeatmap`/`drawMultiHeatmap`：呼叫 `heatmapSVG` 時傳 `clickable:true`；render 後對 svg 掛**委派點擊**：讀 `data-r/data-c`，**依 `state.swap` 反推** `(factorIdx, cancerIdx)`：
  - 不 swap：`factorIdx=+data-r, cancerIdx=+data-c`
  - swap：`factorIdx=+data-c, cancerIdx=+data-r`（因為 swap 後 row=cancer、col=factor）
  - 再 `drawDrillForest(factorIdx, cancerIdx)` 把該癌種 forest 畫進 `#cx-drill`。
- 新 `drawDrillForest(_, cancerIdx)`：用 `cancers[cancerIdx]`，照 `state.model` 重算該癌種的 forest（**可重用** `renderUnivariate`/`renderMultivariate` 的計算，但需改成「畫進指定容器 + 自己的存檔鍵」，或抽出一個 forest 計算函式回傳 items 再各自畫）。建議仿 `correlation.js` 的 `scatterWidget`（自帶控制/存檔、寫進 host 容器）。
- 存檔：heatmap 主鍵存 heatmap；drill 面板自帶 SVG/PNG/CSV（檔名 `cox_forest_{cancer}_{ep}`）。

## 4. 不可破壞的既有行為（硬不變量）
- `heatmapSVG` 對**現有 6 個呼叫端**輸出不變（clickable gated）。寫 smoke：不傳 clickable → 輸出與改前一致。
- **Cox 單癌種 → forest 完全不變**。
- **A=1, B=0 編碼**（HR = advanced/High vs baseline/Low）；`coxPH1`/`coxPH` 引擎**不要動**（已驗證、有對稱性等不變量，見 PROJECT_HANDOFF §7）。
- **SVG 尺寸**：主圖內聯 `style="width:100%;max-width:${W}px;height:auto"`；多 SVG gradient/clip 用唯一 uid。
- **最終圖乾淨**：heatmap 上不寫警告字（淡格/斜線表達）。
- **FDR 每癌種欄內**（別改成整張）。
- swap 轉置與「點擊反推座標」必須一致（上面 §3 的公式）。
- UI/圖/匯出文字**全英文**；註解**繁中**；溝通**繁中**。

## 5. 測試方式
- `node --check js/core/plots.js js/analyses/coxRegression.js`。
- **heatmapSVG smoke**：(1) 傳 `clickable:true` → 含 `class="hm-cell"` + `data-r="0" data-c="1"`；(2) 不傳 → 輸出與改前**完全一致**（比對 viewBox/內容，確認沒影響其他模組）。
- **瀏覽器**（`python -m http.server 8000` → Ctrl+F5）：Cox → 選 2+ 癌種 → uni 與 multi heatmap → **點一格 → 右側出現該癌種 forest**；再點別格 → 更新；heatmap 不消失；heatmap 存檔 vs forest 存檔都對且檔名可分辨；**Swap 後再點 → 對應到正確的 (factor,cancer)**（重點驗 swap 座標反推）。
- 回歸：單癌種 forest 不變；其他模組的 heatmap（survival/immune/correlation/clinicalOverview）外觀不變。

## 6. Part C「完成」的收尾
1. 通過上面測試。
2. **更新** `PROJECT_HANDOFF.md §6⑥`（Cox Regression）補：多癌種 heatmap（uni/multi）+ 點格 drill-down；`§8` 補 `heatmapSVG` 的 clickable opt。
3. **更新 memory** `cox-regression-design`（Part C done）。
4. **刪除本檔 `HANDOFF_PART_C.md`**。
5. commit（建議訊息：`Add click-to-forest drill-down to Cox HR heatmap`）。
