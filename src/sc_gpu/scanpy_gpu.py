import os
import scanpy as sc
try:
    import cuml
    from . import rapids_scanpy_funcs
except ImportError as e:
    print("Please confirm cuml is installed. ")
    raise e


class Globals:

    def __init__(self) -> None:
        self.INITIALIZED = False


class InitializationError(Exception):

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


globals = Globals()


def init(gpu: int = 0):
    # os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu)
    globals.INITIALIZED = True


def _check_init(func):

    def func_with_init_check(*args, **kwargs):
        if not globals.INITIALIZED:
            # raise InitializationError('Please call `scanpy_gpu.init` to initilize first.')
            pass
        return func(*args, **kwargs)

    return func_with_init_check


def _get_rep(adata, use_rep):
    if use_rep is None:
        rep = adata.X
    else:
        rep = adata.obsm[use_rep]
    return rep


@_check_init
def pca(adata, use_rep=None, n_comps=None):
    if n_comps is None:
        n_comps = min(50, adata.shape[0], adata.shape[1])
    pca = cuml.PCA(n_components=n_comps)
    rep = _get_rep(adata, use_rep)
    pca_transformed = pca.fit_transform(rep)
    adata.obsm['X_pca'] = pca_transformed


@_check_init
def umap(adata, use_rep='X_pca', min_dist=0.5, n_neighbors=15):
    umap = cuml.UMAP(n_components=2, min_dist=min_dist, n_neighbors=n_neighbors)
    rep = _get_rep(adata, use_rep)
    umap_2d = umap.fit_transform(rep)
    adata.obsm['X_umap'] = umap_2d


@_check_init
def tsne(adata, use_rep='X_pca', perplexity=30.0):
    tsne = cuml.TSNE(n_components=2, perplexity=perplexity)
    rep = _get_rep(adata, use_rep)
    tsne_2d = tsne.fit_transform(rep)
    adata.obsm['X_tsne'] = tsne_2d


@_check_init
def leiden(adata, resolution=1, use_rep=None, add_key='leiden'):
    if use_rep is None:
        sc.pp.neighbors(adata, use_rep=use_rep)
    adata.obs[add_key] = rapids_scanpy_funcs.leiden(adata, resolution=resolution)