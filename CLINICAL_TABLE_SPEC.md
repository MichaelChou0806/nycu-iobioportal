# 臨床表格式 spec（multivariate / OSCC 用）

> 給離線整理臨床資料用（R/Python，在另一台電腦）。產出一張「寬表」上傳 R2，前端就能用同一套 Cox / dimensions 邏輯。
> **核心原則**：CSV 只放**原始值**；「怎麼切成 baseline/advanced」由 `config/dimensions.*.json` 定義（規則 #3：切法即資料）。

---

## 1. 檔案格式
- 純 CSV、UTF-8、逗號分隔、**一列表頭、一病人一列**。
- 最後 gzip 成 `*.csv.gz` 上 R2（前端用 DecompressionStream 解壓）。
- `manifest.json` 的 `clinical_file` 指到這個 gz、`clinical_fields` 列出所有欄名（你現有 pipeline 已這樣做）。

## 2. 必要欄
| 欄 | 規則 |
|---|---|
| `patient_id` | join key，**必須和 `samples.csv` 的 patient_id 完全一致**。沒有對應 tumor sample 的列不會被用到。 |

## 3. 存活 endpoint（Cox 需要）
每個 endpoint 兩欄：`<EP>` = 事件指標（1=event, 0=censored）、`<EP>.time` = 追蹤時間（**天**）。
- 至少 `OS` + `OS.time`。
- 有就放 `DSS`/`DSS.time`、`DFI`/`DFI.time`、`PFI`/`PFI.time` —— 前端自動偵測哪些可用，下拉只顯示有資料的。

## 4. 臨床因子
- 一因子一欄，**原始值，不要 dummy、不要先分箱**。
  - 類別：保留原字串（`T2a`、`Stage IVA`、`POSITIVE`、`G3`、`YES`/`NO`），**同一欄拼法一致**。
  - 連續：放數字（`age=63`、`pack_years=40`）。
- **缺值 = 留空白**（前端把空白當 missing → 該病人在該因子被排除；multivariate 進不了 complete-case）。

## 5. 其他欄
- `redaction`：有被撤回的病人就保留此欄（非空者會被前端 `patientsInScope` 排除）。
- 一致性：patient_id 格式一致、無重複列、同欄不混單位。

## 6. 因子怎麼變成「可選」（重要）
一個 CSV 欄位**不會自動變成 Cox 可選因子**。要讓它可選、可切：在一份 dimensions JSON 加一筆（格式同 `config/dimensions.tcga.json`）：
- `field` 指向 CSV 欄名。
- `type: "binary"`（天然兩組，給 `baseline`/`advanced` 的值清單或 `numericSplit` 的 cutoff）或 `type: "ordinal"`（多級，每級指派 baseline/advanced/ignore）。
- OSCC 另開一份（如 `config/dimensions.oscc.json`），在 `datasets.json` 的該筆用 `dimensionsUrl` 指過去；相同維度沿用、TCGA 沒有的（如 ENE）在 OSCC 版加上。

→ **CSV 給原始欄；dimensions JSON 定義怎麼切。** 任意原始欄的「深撈」是之後 factor picker 的事。

## 7. 範例（空白 = 缺值）
```csv
patient_id,OS,OS.time,DSS,DSS.time,gender,age,pathologic_T,node_status,pathologic_stage,perineural_invasion_present,hpv_status
TCGA-XX-1234,1,612,1,612,MALE,63,T2a,N+,Stage IVA,YES,Positive
TCGA-XX-5678,0,1430,0,1430,FEMALE,71,T1,N-,Stage I,NO,Negative
TCGA-XX-9012,1,205,,,MALE,,T4,N+,,,
```

## 8. OSCC 接入注意
- 免疫：OSCC 跑自己的 deconvolution，`build_immune.py` 設 `KEY_MODE="exact"`（OSCC sample_id 非 TCGA barcode，直接對）。
- 補 ENE 維度（OSCC 專屬）。
- 地雷：batch effect（最大風險，記錄 batch、相同 pipeline、考慮 ComBat）；小樣本 power（~100 tumor+30 normal，quartile 後每組 ~25，Cox/KM 不穩）。

## 9. git / 資料處理
- **真實臨床表是資料，絕對不要 commit 進這個網站 repo**（它會上 Cloudflare Pages＝公開，且是病人資料）。整理在 pipeline 那邊、gitignore data、走 R2。
- 建議在各 repo 的 `.gitignore` 加 `*.csv` / `*.csv.gz` 當保險。
- 本檔（spec）是文件，可以 commit。
