# sc_gpu
single cell utilities on gpu

Accelerates scanpy functions `sc.tl.umap`, `sc.tl.tsne`, `sc.tl.leiden` and scib functions `scib.metrics.silhouette`, `scib.metrics.silhouette_batch` by 3-100 folds. 

see [examples](experiments/examples.ipynb).

requires cuml environment; see [rapids](https://rapids.ai/start.html#get-rapids).