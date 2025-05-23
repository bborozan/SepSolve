import numpy as np
import pandas

from .sepsolve_base import MarkerGeneLPSolver
from .sepsolve_fixed import SepSolveFixed

def __get_markers__internal(data, labels, num_markers, s=0.4, ilp=False):
    base = MarkerGeneLPSolver(data, labels, num_markers, ilp=ilp)
    solver = SepSolveFixed(base, s)
    
    x, betas, obj = base.Solve(solver)
    
    return base.ranking(x)

def __process_labels(labels):
    # convert pandas series into numpy array
    if isinstance(labels, pandas.core.series.Series):
        return labels.astype("category").cat.codes.to_numpy()
    else:
        categorical = pandas.Categorical(labels)
        return np.array(categorical.codes)
        
def get_markers(data, labels, num_markers, c=0.4, ilp=False):
    """Return indices of `num_markers` marker genes providing *c*-separation.

    Parameters
    ----------
    data        : ndarray, shape (cells, genes)
        Pre‑processed expression matrix.
    labels      : array‑like, shape (cells,)
        Integer or string cell‑type annotations.
    num_markers : int
        Number of markers to select.
    c      : float, default 0.4
        Separation strictness. Larger → stricter.
    ilp : bool, default False
        Use ILP solver if True; otherwise, use LP relaxation (faster, approximate).

    Returns
    -------
    list of int
        Indices of selected marker genes.
    """

    if data.size == 0:
        raise ValueError(f"Data size is zero.")
    
    # convert labels to integers
    lab = __process_labels(labels)

    # check if there is a label for every cell
    if len(lab) != data.shape[0]:
        # check if data has to be transposed
        # we want cells as rows and genes as columns
        raise ValueError(f"Label vector has length {len(lab)}, while the data matrix is of shape {data.shape}.")
    
    return __get_markers__internal(data, lab, num_markers, s=c, ilp=ilp)
