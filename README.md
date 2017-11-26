# Gas_adsorption_lib
Calculate BJH pore size distribution

## How to install 
simple example 
```python
from adsorp_libary.BJH_calculation import BJH_method
solution = BJH_method(gas_type=gas_type)
solution.fit_transform(p=p_rels, q=q_abs)

```