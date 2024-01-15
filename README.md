
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stemv

<!-- badges: start -->
<!-- badges: end -->

このパッケージは単木の胸高直径、樹高および樹種と地域から立木幹材積を計算する関数を提供します。  
計算方法は林野庁（1970）「立木幹材積表 東日本編」「立木幹材積表
西日本編」の書籍とその修正を 提案した細田ら（2010）森林計画学会誌
44:23-39 の論文に基づいています。  
上記の細田ら(2010)
の論文に従った計算を行うExcelのユーザ定義関数がすでに配布されており、
このパッケージはそれに類似しています（ただし、下記に記す通り差異があります）。

## Installation

`devtools`を利用して最新版の`stemv`パッケージをGitHubからインストールできます。

``` r
devtools::install_github("dulvrq/stemv")
library(stemv)
```

## Example

利用する立木幹材積式の種類（地域と樹種）がわかっている場合、胸高直径（DBH）と樹高から材積を計算できます。

``` r
library(stemv)

Name <- "東京スギ"
D <- 15.1
H <- 20.2
stemVolume(Name, D, H)
#> [1] 0.1893091
```

地域と樹種から立木幹材積式の種類を特定できます。地域名は都道府県支庁名または旧営林局です。
北海道については、森林生態系多様性基礎調査データを処理できるよう森林計画区も指定できるように設定しています。

``` r
Region <- "高知県"
Spp <- "ヒノキ"
volumeName(Region, Spp)
#> [1] "高知ヒノキ"
```

`volumeName()`を利用する場合、樹種はカタカナで記入してください。ひらがな等にすると全て広葉樹になります。

``` r
## correct
volumeName("福岡", "カラマツ")
#> [1] "名古屋カラマツ"

## incorrect
volumeName("福岡", "からまつ")
#> [1] "熊本広葉樹2類"
```

樹種のリストは森林総研が配布している「幹材積計算プログラム」を参照してください。
リストにないものは全て広葉樹になります。

``` r
## Banana is not in the list
volumeName("東京", "バナナ")
#> [1] "東京広葉樹"
```

リストにない地域を指定した場合、エラーになります。
ただし、特定の文字列を返すようにすることはできます。

``` r
## cause error because America is not in the list
volumeName("アメリカ", "トドマツ")
#> Error in volumeNameSingle(Region, Spp, RS): No such region.

## return specific string if the region does not much ones in the list
volumeName("アメリカ", "トドマツ", name_invalid = "合致しない")
#> Error in volumeName("アメリカ", "トドマツ", name_invalid = "合致しない"): unused argument (name_invalid = "合致しない")
```

立木幹材積の計算はベクトルでも実行できます。

``` r
dt <- data.frame(Region = c("宗谷支庁", "愛媛県", "東京", "青森県", "熊本"),
                 Spp = c("トドマツ", "ヒノキ", "スギ", "カラマツ", "ヒノキ"),
                 D = c(21.8, 55.1, 31.8, 43.4, 39.0),
                 H = c(11.4, 22.7, 23.4, 19.5, 26.7))
dt |> 
  dplyr::mutate(Name = volumeName(Region, Spp),
                V = stemVolume(Name, D, H))
#>     Region      Spp    D    H         Name         V
#> 1 宗谷支庁 トドマツ 21.8 11.4 旭川トドマツ 0.2287186
#> 2   愛媛県   ヒノキ 55.1 22.7   高知ヒノキ 2.2230540
#> 3     東京     スギ 31.8 23.4     東京スギ 0.8597140
#> 4   青森県 カラマツ 43.4 19.5 秋田カラマツ 1.3602628
#> 5     熊本   ヒノキ 39.0 26.7   熊本ヒノキ 1.4668625
```

幹材積の計算では、細田ら(2010)の方法に従い必要な箇所は材積式の各接合部に移動平均を適用しています。
もし、移動平均を適用しない場合の計算結果を求めたい場合、`stemVolumeSingle`関数を利用してください。
ただし、この関数はベクトルを処理できないので単一の値のみ渡してください。
複数の立木を処理したい場合はfor構文か`purrr::map`系列の関数を利用してください

``` r
## apply moving average adjustment
stemVolumeSingle("函館エゾマツ", 51, 25, off_adj = FALSE)
#> [1] 2.29211
## do not apply moving average adjustment
stemVolumeSingle("函館エゾマツ", 51, 25, off_adj = TRUE)
#> [1] 2.35948
```

``` r
## if multiple values are passed (very slow for large amount of data)
purrr::pmap(list(Name = volumeName(dt$Region, dt$Spp),
                 D = dt$D,
                 H = dt$H), stemVolumeSingle, off_adj = TRUE) |> unlist()
#> [1] 0.2200923 2.2230540 0.8597140 1.3602628 1.4668625
```

## Differnces from calculation program (Excel version)

「幹材積計算プログラム」とこのパッケージは独立しており、いくつかの点で差異があります。

- 幹材積計算プログラムにおける「札幌トドマツ」、「高知天然スギ」、「青森アカマツ」に対する
  追加の修正を反映していない。
- 「青森広葉樹」、「高知広葉樹」の5点移動平均の当てはめ範囲が異なる。
- 胸高形数法における四捨五入計算の際、一部で浮動小数点の誤差が発生する（北海道針葉樹）。
- 一部の計算における係数の有効数字が異なる。

## Reference

細田和男, 光田 靖, 家原敏郎 (2010)
現行立木幹材積表と材積式による計算値との相違およびその修正方法」
森林計画学会誌 44: 23-39. <https://doi.org/10.20659/jjfp.44.2_23>

森林総合研究所「幹材積計算プログラム」
<https://www.ffpri.affrc.go.jp/database/stemvolume/index.html>
