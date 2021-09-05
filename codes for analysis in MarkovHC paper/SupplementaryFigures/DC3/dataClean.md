```R
setwd('/data02/zywang/MarkovHC/supplementaryFigures/DC3(scATACSeq+scRNASeq)/')
```


```R
matrix <- read.table(file = './GSE107651_scATAC-seq_RA_D4.txt', sep='\t', header = T)
```


```R
matrix$chromosome <- as.character(matrix$chromosome)
```


```R
head(matrix)
```


<table>
<caption>A data.frame: 6 × 99</caption>
<thead>
	<tr><th></th><th scope=col>chromosome</th><th scope=col>start</th><th scope=col>end</th><th scope=col>scATAC.seq_RA_D4_S1</th><th scope=col>scATAC.seq_RA_D4_S2</th><th scope=col>scATAC.seq_RA_D4_S3</th><th scope=col>scATAC.seq_RA_D4_S4</th><th scope=col>scATAC.seq_RA_D4_S5</th><th scope=col>scATAC.seq_RA_D4_S6</th><th scope=col>scATAC.seq_RA_D4_S7</th><th scope=col>⋯</th><th scope=col>scATAC.seq_RA_D4_S87</th><th scope=col>scATAC.seq_RA_D4_S88</th><th scope=col>scATAC.seq_RA_D4_S89</th><th scope=col>scATAC.seq_RA_D4_S90</th><th scope=col>scATAC.seq_RA_D4_S91</th><th scope=col>scATAC.seq_RA_D4_S92</th><th scope=col>scATAC.seq_RA_D4_S93</th><th scope=col>scATAC.seq_RA_D4_S94</th><th scope=col>scATAC.seq_RA_D4_S95</th><th scope=col>scATAC.seq_RA_D4_S96</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>⋯</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1</td><td>3107081</td><td>3108222</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>2</th><td>chr1</td><td>3110845</td><td>3111505</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>3</th><td>chr1</td><td>3345258</td><td>3345557</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>chr1</td><td>3452492</td><td>3452796</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>5</th><td>chr1</td><td>3569409</td><td>3570892</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>chr1</td><td>3931685</td><td>3931917</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
for(i in 1:nrow(matrix)){
    rownames(matrix)[i] <- paste(matrix[i,1:3], collapse = "_")
}
```


```R
dim(matrix)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>105095</li><li>99</li></ol>




```R
head(matrix)
```


<table>
<caption>A data.frame: 6 × 99</caption>
<thead>
	<tr><th></th><th scope=col>chromosome</th><th scope=col>start</th><th scope=col>end</th><th scope=col>scATAC.seq_RA_D4_S1</th><th scope=col>scATAC.seq_RA_D4_S2</th><th scope=col>scATAC.seq_RA_D4_S3</th><th scope=col>scATAC.seq_RA_D4_S4</th><th scope=col>scATAC.seq_RA_D4_S5</th><th scope=col>scATAC.seq_RA_D4_S6</th><th scope=col>scATAC.seq_RA_D4_S7</th><th scope=col>⋯</th><th scope=col>scATAC.seq_RA_D4_S87</th><th scope=col>scATAC.seq_RA_D4_S88</th><th scope=col>scATAC.seq_RA_D4_S89</th><th scope=col>scATAC.seq_RA_D4_S90</th><th scope=col>scATAC.seq_RA_D4_S91</th><th scope=col>scATAC.seq_RA_D4_S92</th><th scope=col>scATAC.seq_RA_D4_S93</th><th scope=col>scATAC.seq_RA_D4_S94</th><th scope=col>scATAC.seq_RA_D4_S95</th><th scope=col>scATAC.seq_RA_D4_S96</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>⋯</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>chr1_3107081_3108222</th><td>chr1</td><td>3107081</td><td>3108222</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3110845_3111505</th><td>chr1</td><td>3110845</td><td>3111505</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3345258_3345557</th><td>chr1</td><td>3345258</td><td>3345557</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3452492_3452796</th><td>chr1</td><td>3452492</td><td>3452796</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3569409_3570892</th><td>chr1</td><td>3569409</td><td>3570892</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3931685_3931917</th><td>chr1</td><td>3931685</td><td>3931917</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
matrix <- matrix[,4:ncol(matrix)]
```


```R
head(matrix)
```


<table>
<caption>A data.frame: 6 × 96</caption>
<thead>
	<tr><th></th><th scope=col>scATAC.seq_RA_D4_S1</th><th scope=col>scATAC.seq_RA_D4_S2</th><th scope=col>scATAC.seq_RA_D4_S3</th><th scope=col>scATAC.seq_RA_D4_S4</th><th scope=col>scATAC.seq_RA_D4_S5</th><th scope=col>scATAC.seq_RA_D4_S6</th><th scope=col>scATAC.seq_RA_D4_S7</th><th scope=col>scATAC.seq_RA_D4_S8</th><th scope=col>scATAC.seq_RA_D4_S9</th><th scope=col>scATAC.seq_RA_D4_S10</th><th scope=col>⋯</th><th scope=col>scATAC.seq_RA_D4_S87</th><th scope=col>scATAC.seq_RA_D4_S88</th><th scope=col>scATAC.seq_RA_D4_S89</th><th scope=col>scATAC.seq_RA_D4_S90</th><th scope=col>scATAC.seq_RA_D4_S91</th><th scope=col>scATAC.seq_RA_D4_S92</th><th scope=col>scATAC.seq_RA_D4_S93</th><th scope=col>scATAC.seq_RA_D4_S94</th><th scope=col>scATAC.seq_RA_D4_S95</th><th scope=col>scATAC.seq_RA_D4_S96</th></tr>
	<tr><th></th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>⋯</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>chr1_3107081_3108222</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3110845_3111505</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3345258_3345557</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3452492_3452796</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3569409_3570892</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>chr1_3931685_3931917</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>⋯</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
write.table(matrix, file='./GSE107651_scATAC-seq_RA_D4_clean.txt', quote = FALSE, row.names = T, col.names = T)
```


```R
RNAmatrix <- read.table(file = './GSE115968_scRNA-seq_RA_D4.txt', header = T)
```


```R
dim(RNAmatrix)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>21972</li><li>465</li></ol>




```R
group <- RNAmatrix[, 1] 
RNAmatrix <- aggregate(RNAmatrix[, 2:ncol(RNAmatrix)], by = list(group), FUN = sum)
```


```R
dim(RNAmatrix)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>21970</li><li>465</li></ol>




```R
rownames(RNAmatrix) <- RNAmatrix$Symbol
```


```R
head(RNAmatrix)
```


<table>
<caption>A data.frame: 6 × 465</caption>
<thead>
	<tr><th></th><th scope=col>Group.1</th><th scope=col>WW31.1.S100_S190_L007</th><th scope=col>WW31.1.S101_S201_L007</th><th scope=col>WW31.1.S102_S202_L007</th><th scope=col>WW31.1.S103_S203_L007</th><th scope=col>WW31.1.S104_S204_L007</th><th scope=col>WW31.1.S105_S205_L007</th><th scope=col>WW31.1.S106_S206_L007</th><th scope=col>WW31.1.S107_S207_L007</th><th scope=col>WW31.1.S109_S209_L007</th><th scope=col>⋯</th><th scope=col>WW31.2.S90_S180_L008</th><th scope=col>WW31.2.S91_S191_L008</th><th scope=col>WW31.2.S92_S192_L008</th><th scope=col>WW31.2.S94_S194_L008</th><th scope=col>WW31.2.S95_S195_L008</th><th scope=col>WW31.2.S96_S196_L008</th><th scope=col>WW31.2.S97_S197_L008</th><th scope=col>WW31.2.S98_S198_L008</th><th scope=col>WW31.2.S99_S199_L008</th><th scope=col>WW31.2.S9_S19_L008</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>0610009O20Rik</td><td>  0.48</td><td>  9.25</td><td>  5.43</td><td>  0.52</td><td>  0.00</td><td>  0.50</td><td>  0.00</td><td>  1.58</td><td>  2.33</td><td>⋯</td><td> 11.03</td><td>  0.00</td><td>  0.00</td><td>  0.00</td><td>  0.00</td><td>  0.00</td><td> 0.00</td><td>  3.53</td><td>  0.00</td><td>  0.00</td></tr>
	<tr><th scope=row>2</th><td>0610010F05Rik</td><td> 31.35</td><td> 13.02</td><td> 41.43</td><td>  0.00</td><td>  0.00</td><td>  2.04</td><td>  0.00</td><td> 26.75</td><td> 57.05</td><td>⋯</td><td> 37.53</td><td>  0.00</td><td> 28.47</td><td>  0.00</td><td> 71.88</td><td> 19.90</td><td>76.98</td><td>  0.75</td><td> 17.16</td><td>  0.00</td></tr>
	<tr><th scope=row>3</th><td>0610010K14Rik</td><td> 95.75</td><td>226.57</td><td>117.39</td><td>253.78</td><td>148.96</td><td> 96.35</td><td>166.16</td><td>142.22</td><td>190.83</td><td>⋯</td><td>206.30</td><td>209.23</td><td>113.25</td><td>244.98</td><td>223.88</td><td>112.98</td><td>77.42</td><td>293.62</td><td>180.65</td><td>272.36</td></tr>
	<tr><th scope=row>4</th><td>0610012G03Rik</td><td>141.58</td><td> 34.29</td><td> 39.30</td><td>155.83</td><td>106.80</td><td>  0.00</td><td> 74.32</td><td>147.53</td><td>104.36</td><td>⋯</td><td> 71.24</td><td>116.11</td><td> 27.85</td><td> 74.46</td><td>105.15</td><td> 73.96</td><td> 0.00</td><td> 47.76</td><td> 44.15</td><td> 44.43</td></tr>
	<tr><th scope=row>5</th><td>0610030E20Rik</td><td>  0.00</td><td>  0.00</td><td> 19.84</td><td>  0.31</td><td>  0.00</td><td>  0.00</td><td>  8.08</td><td>  0.00</td><td>  0.00</td><td>⋯</td><td>  0.00</td><td> 51.97</td><td>  0.62</td><td>  0.00</td><td>  0.00</td><td> 55.76</td><td>63.57</td><td> 15.51</td><td>  0.00</td><td> 19.59</td></tr>
	<tr><th scope=row>6</th><td>0610037L13Rik</td><td>348.92</td><td>222.96</td><td> 88.13</td><td> 94.59</td><td>168.82</td><td>279.19</td><td>209.39</td><td>111.89</td><td>198.31</td><td>⋯</td><td> 89.84</td><td>247.13</td><td> 91.64</td><td> 62.02</td><td> 95.84</td><td>230.68</td><td>65.63</td><td>181.89</td><td> 86.70</td><td>176.14</td></tr>
</tbody>
</table>




```R
rownames(RNAmatrix) <- RNAmatrix[,1]
```


```R
RNAmatrix <- RNAmatrix[,2:ncol(RNAmatrix)]
```


```R
write.table(RNAmatrix, file='./GSE115968_scRNA-seq_RA_D4_clean.txt', quote = FALSE, row.names = T, col.names = T)
```


```R

```
