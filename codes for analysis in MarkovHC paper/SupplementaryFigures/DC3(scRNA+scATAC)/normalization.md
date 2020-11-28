```python
import numpy as np
import pandas as pd
```


```python
def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df
```


```python
E = pd.read_csv('/data02/zywang/MarkovHC/DC3/GSE107651_scATAC-seq_RA_D4_clean.txt',sep=' ',header=0)
E = pd.DataFrame(E)
```


```python
E = quantileNormalize(E) 
```


```python
RNA = pd.read_csv('/data02/zywang/MarkovHC/DC3/GSE115968_scRNA-seq_RA_D4_clean.txt',sep=' ',header=0)
RNA = pd.DataFrame(RNA)
```


```python
RNA = quantileNormalize(RNA) 
```


```python
RNA.to_csv("/data02/zywang/MarkovHC/DC3/scRNA.txt", sep=" ")
E.to_csv("/data02/zywang/MarkovHC/DC3/scATAC.txt", sep=" ")
```


```python

```
