# DM_final

## 執行環境

以下為 2026-11-10 從本機 R 環境查得的版本，記錄 SCOTH 重現程式使用的直接依賴。

- **R**：4.6.2（2026-10-29 ucrt）
- **平台**：Windows 64 位元（`x86_64-w64-mingw32`）

| R 套件 | 已安裝版本 |
| --- | --- |
| scoth | 1.0 |
| dplyr | 1.2.1 |
| ggplot2 | 4.0.3 |
| grid | 4.6.2 |
| gridExtra | 2.3.1 |
| here | 1.0.2 |
| jsonlite | 2.0.0 |
| logger | 0.4.3 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| testthat | 3.3.2 |
| pkgload | 1.5.3 |

`grid` 隨 R 提供；`scoth` 為本機安裝的原作者套件。`setup/install_packages.R` 未鎖定套件版本，此表供重現環境參考。

SOTA 比較腳本 `project0_reproduction/experiments/main_sota.R` 所需的額外 R 套件目前尚未安裝：`breakfast`、`changepoint`、`changepoint.np`、`ecp`、`genlasso`、`ocp`、`reticulate`、`robseg`、`scmamp`。

# 專案目錄
```
final/
├── setup/
│   └── install_packages.R       # 安裝 R 套件
├── project0_reproduction/       # 原論文重現
│   ├── R/                       # 演算法與評估程式
│   ├── experiments/             # 實驗執行腳本
│   ├── public-data/             # 資料集與標註
│   ├── tests/                  ㄒㄩ # 單元測試
│   ├── man/                     # R 函式說明文件
│   ├── output-data/             # 實驗中間檔與結果
│   ├── results/                 # 整理後的圖表
│   ├── logs/                    # 執行紀錄
│   └── README_original_authors.md # 原作者說明
├── reference/
│   └── scoth-segmentation/      # 原作者程式碼參考副本
```
