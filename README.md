# Gas_adsorption_lib
Calculate BJH pore size distribution

## How to install 
simple example 
```python

# import pakages
import pandas as pd
from adsorp_libary.BJH_calculation import BJH_method

# Get data
filename = 'adsorp_libary/data/test_isos.xlsx'
sheetname = '3_53_N2'
df_iso = pd.read_excel(filename,sheetname) 
p_rels = df_iso['p_rels'].values
q_abs =  df_iso['q_abs'].values

p_rels = df_iso['p_rels'].values
q_abs =  df_iso['q_abs'].values

# Call BJH class 
solution = BJH_method(gas_type='N2')
solution.fit_transform(p=p_rels, q=q_abs)

# Results
Davg = solution.Davg # pore diameter 
Vp = solution.Vp # volume distribution
Davg, Vp

```
For details, please check testing jupyter notebook