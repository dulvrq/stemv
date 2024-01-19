
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stemv

<!-- badges: start -->
<!-- badges: end -->

このパッケージは単木の胸高直径、樹高および樹種と地域（日本国内）から立木幹材積を計算する関数を提供します。
計算方法は林野庁（1970）「立木幹材積表 東日本編」「立木幹材積表
西日本編」の書籍とその修正を提案した細田ら（2010）森林計画学会誌
44:23-39 の論文に基づいています。  
上記の細田ら(2010)
の論文に従った計算を行うExcelのユーザ定義関数がすでに配布されており、
このパッケージはそれに類似しています（ただし、下記に示す通り差異があります）。

## Installation

`devtools`を利用して最新版の`stemv`パッケージをGitHubからインストールできます。

``` r
devtools::install_github("dulvrq/stemv")
library(stemv)
```

## Example

利用する立木幹材積式の種類（地域と樹種の組み合わせ）がわかっている場合、胸高直径（DBH）と樹高（H）から材積を計算できます。

``` r
library(stemv)

Name <- "東京スギ"
D <- 15.1
H <- 20.2
stemVolume(Name, D, H)
#> [1] 0.1893091
```

地域と樹種から立木幹材積式の種類を特定できます。地域名は都道府県支庁名または旧営林局名です。
北海道については、森林生態系多様性基礎調査データの処理のため森林計画区名も指定できるように設定しています。

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
```

``` r
## return specific string if the region does not much any one in the list
volumeNameSingle("アメリカ", "トドマツ", name_invalid = "合致しない")
#> [1] "合致しない"
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
もし、移動平均を適用しない計算結果を求めたい場合、`stemVolumeSingle`関数を利用してください。
ただし、この関数はベクトルを処理できないので単一の値のみ渡してください。
複数の立木を処理したい場合はfor構文か`purrr::map`系列の関数を利用してください。

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

## Differences from the calculation program (Excel version)

「幹材積計算プログラム」（Excel版）とこのパッケージは独立しており、いくつかの点で差異があります。

1.  幹材積計算プログラムにおける「札幌トドマツ」、「高知天然スギ」、「青森アカマツ」に対する追加の修正を反映していない。
2.  「青森広葉樹」、「高知広葉樹」の5点移動平均の当てはめ範囲が異なる。
3.  一部の計算における係数の有効数字が異なる。

各材積計算式におけるDBH: 0-200 cm (1 cmごと) および H: 0-50 m (1 mごと)
の各組み合わせについて、幹材積計算プログラムと
`stemVolume()`との材積の計算結果を比較すると以下のようになります。

- 誤差が0.0001m<sup>3</sup>以上になる組み合わせは、全体の1.4%で0.01m<sup>3</sup>以上になる組み合わせは0.2%です。
- このうち、0.01m<sup>3</sup>以上の誤差は全て上記1. ~ 2.
  に起因し、0.0001m<sup>3</sup>以上の誤差は全て上記1. ~ 3.に起因します。

<table>
<caption>
表1. 差が存在する材積式における誤差
</caption>
<thead>
<tr>
<th style="text-align:left;">
計算式名
</th>
<th style="text-align:right;">
組み合わせ数
</th>
<th style="text-align:right;">
差≥0.01の数
</th>
<th style="text-align:right;">
差≥0.0001の数
</th>
<th style="text-align:right;">
平均絶対誤差
</th>
<th style="text-align:right;">
最大誤差
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
札幌トドマツ
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
428
</td>
<td style="text-align:right;">
459
</td>
<td style="text-align:right;">
0.0033
</td>
<td style="text-align:right;">
0.1438
</td>
</tr>
<tr>
<td style="text-align:left;">
函館エゾマツ
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
6925
</td>
<td style="text-align:right;">
0.0005
</td>
<td style="text-align:right;">
0.0009
</td>
</tr>
<tr>
<td style="text-align:left;">
青森アカマツ
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
656
</td>
<td style="text-align:right;">
850
</td>
<td style="text-align:right;">
0.0038
</td>
<td style="text-align:right;">
0.2963
</td>
</tr>
<tr>
<td style="text-align:left;">
青森広葉樹
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0.0002
</td>
<td style="text-align:right;">
0.0435
</td>
</tr>
<tr>
<td style="text-align:left;">
長野カラマツ
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3336
</td>
<td style="text-align:right;">
0.0001
</td>
<td style="text-align:right;">
0.0005
</td>
</tr>
<tr>
<td style="text-align:left;">
高知天然スギ
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
353
</td>
<td style="text-align:right;">
416
</td>
<td style="text-align:right;">
0.0028
</td>
<td style="text-align:right;">
0.6360
</td>
</tr>
<tr>
<td style="text-align:left;">
高知広葉樹
</td>
<td style="text-align:right;">
10251
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
0.0002
</td>
<td style="text-align:right;">
0.0556
</td>
</tr>
</tbody>
</table>

## Differences from the stem volume table

「立木幹材積表
東日本編・西日本編」における幹材積表と比較すると以下の通りになります。

- 小数点第2位を四捨五入して材積が一致しない組み合わせは全体の4.7%です。
- 誤差が0.01m<sup>3</sup>以上になる組み合わせは全体の1.2%です。

<table>
<caption>
表2. 一致しない数が多い10種類の材積式における誤差
</caption>
<thead>
<tr>
<th style="text-align:left;">
計算式名
</th>
<th style="text-align:right;">
組み合わせ数
</th>
<th style="text-align:right;">
一致しない数
</th>
<th style="text-align:right;">
差≥0.01の数
</th>
<th style="text-align:right;">
平均絶対誤差
</th>
<th style="text-align:right;">
最大誤差
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
青森ヒバ
</td>
<td style="text-align:right;">
1313
</td>
<td style="text-align:right;">
408
</td>
<td style="text-align:right;">
107
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
札幌トドマツ
</td>
<td style="text-align:right;">
1714
</td>
<td style="text-align:right;">
303
</td>
<td style="text-align:right;">
134
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
北見エゾマツ
</td>
<td style="text-align:right;">
1206
</td>
<td style="text-align:right;">
244
</td>
<td style="text-align:right;">
112
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
青森アカマツ
</td>
<td style="text-align:right;">
803
</td>
<td style="text-align:right;">
222
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
高知天然スギ
</td>
<td style="text-align:right;">
2251
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
172
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.31
</td>
</tr>
<tr>
<td style="text-align:left;">
青森広葉樹
</td>
<td style="text-align:right;">
1972
</td>
<td style="text-align:right;">
198
</td>
<td style="text-align:right;">
78
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
札幌エゾマツ
</td>
<td style="text-align:right;">
1696
</td>
<td style="text-align:right;">
183
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
関東中部ツガ
</td>
<td style="text-align:right;">
1592
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
函館ヒバ
</td>
<td style="text-align:right;">
1164
</td>
<td style="text-align:right;">
179
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
高知天然ヒノキ
</td>
<td style="text-align:right;">
1985
</td>
<td style="text-align:right;">
177
</td>
<td style="text-align:right;">
103
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
</tbody>
</table>

## Reference

細田和男, 光田 靖, 家原敏郎 (2010)
現行立木幹材積表と材積式による計算値との相違およびその修正方法」
森林計画学会誌 44: 23-39. <https://doi.org/10.20659/jjfp.44.2_23>

森林総合研究所「幹材積計算プログラム」
<https://www.ffpri.affrc.go.jp/database/stemvolume/index.html>
