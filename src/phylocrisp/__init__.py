from ._consensus import Tree_with_support, TreeList_with_support, transfer_support, unnormalized_transfer_support
from ._simulate_tree import birthdeath_sampling, normalize_tree_with_lognormal
from ._greedy import STDGreedyConsensus, SQDGreedyConsensus, SUTDGreedyConsensus
from ._metrics import SUTD_fp_fn, STD_fp_fn, SQD_fp_fn, SBD_fp_fn
from ._booster import load_booster
from ._c_pruning import c_prune

_VISUALIZE_EXPORTS = (
    "plot_example_func",
    "get_support",
    "get_support_from_NeXML",
    "write_consensus_NeXML",
    "read_consensus_NeXML",
)


def __getattr__(name):
    """Lazily import visualization helpers so ete4 is only required when used."""
    if name in _VISUALIZE_EXPORTS:
        try:
            from . import _visualize as _viz
        except ModuleNotFoundError as e:
            if e.name and e.name.startswith("ete4"):
                raise ModuleNotFoundError(
                    "Visualization helpers require the optional dependency 'ete4'. "
                    "Install it with `pip install \"PhyloCRISP[viz]\"`."
                ) from e
            raise
        value = getattr(_viz, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [s for s in dir() if not s.startswith('_')]
__all__.extend(_VISUALIZE_EXPORTS)
