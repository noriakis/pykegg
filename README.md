# pykegg
Analyze and visualize KEGG information using network approach.

Using `biopython`, `igraph` and [`plotnine`](https://github.com/has2k1/plotnine), parse the KGML into `igraph` and easily plot the relevant information contained in KEGG Pathway in Python environment.

```python
import requests_cache
import numpy as np
from PIL import Image
import pykegg

## Be sure to cache all the downloaded files to avoid recursive querying
requests_cache.install_cache('pykegg_cache')
g = pykegg.KGML_graph(pid="hsa03460")

## Export as igraph
gra = g.get_graph()
print(gra)

## Overlay to raw image
nds = g.get_nodes()
nds["lfc"] = np.arange(-2,2,0.01)[0:nds.shape[0]]
nds["color"] = nds.lfc.apply(lambda x: round(x,3)).replace(pykegg.color_grad(minimum=min(nds.lfc), maximum=max(nds.lfc),seq=0.01))
Image.fromarray(pykegg.overlay_opencv_image(nds, pid="hsa03460"))

## Plot using plotnine
# options.figure_size = (7,5)
pykegg.plot_kegg_pathway_plotnine(g, node_x_nudge=25, label_size=5, show_label="gene")
```

## [Documentation](https://pykegg.readthedocs.io/)

## TODO
- [ ] The function for converting identifiers (especially for ENTREZID <=> SYMBOL) without connection to servers, like using [`genomic-features`](https://genomic-features.readthedocs.io/en/latest/).
- [ ] Parallel edge support by nudging x and y position