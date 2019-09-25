## RNA velocity analysis

### Data
This ipython notebook requires `velocyto` prepared `loom` files for TiO2 and Asbestos
datasets:
  * `../../SC_14.loom`: TiO2 dataset
  * `../../SC_15.loom`: Asbestos

We prepared them with the following command using Mouse genome reference version GRCm38.p5

```
velocyto run10x <path-to-cellranger-output> <genome-gtf file>
```

It also requires `csv` files with coordinates and identities of cells, which need to be exported
from Seurat objects prepared with this script: https://github.com/NUPulmonary/JoshiWatanabe2019/blob/master/Joshi_Watanabe_ERJ_2019_asbestos/Joshi_Watanabe_ERJ_2019_asbestos.R

  * `../../macrophages-*.csv` from `Macrophages2` object (`Macrophages_clean.Robj` file)
  * `../../at2-*.csv` from `AT2_clean` objects (`AT2_clean.Robj` file)

### Dependencies
We used Python version 3.6.9

Dependencies are fixed in `Pipfile.lock`, please use `pipenv` to install them

