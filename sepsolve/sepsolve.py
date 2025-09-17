import numpy as np
import pandas

from .sepsolve_base import MarkerGeneLPSolver
from .sepsolve_fixed import SepSolveFixed
from .sepsolve_param_opt import __optimise_c_internal

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
    """Selects marker genes from the provided gene expression data based on the specified separation parameter c and the number of markers to select.

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

def optimize_c(data, labels, num_markers, start=0.2, end=1.0, step_size=0.025, verbose=False):
    """
    Optimize the separation parameter `c` for marker gene selection.

    This function evaluates a range of candidate `c` values to determine
    which one provides the best separation between clusters. For each
    candidate `c`, the function selects `num_markers` marker genes 
    using the `get_markers` method and evaluates cluster separation.

    Parameters
    ----------
    data : array-like, shape (n_cells, n_genes)
        Preprocessed gene expression matrix.
    labels : array-like, shape (n_cells,)
        Cluster or cell type annotations for each cell.
    num_markers : int
        Number of marker genes to select for each candidate `c`.
    start : float, optional, default=0.2
        Starting value of `c` to evaluate.
    end : float, optional, default=1.0
        Ending value of `c` to evaluate.
    step_size : float, optional, default=0.025
        Step size for scanning `c` values between `start` and `end`.
    verbose : bool, optional, default=False
        If True, prints progress and evaluation results for each `c`.

    Returns
    -------
    optimal_c : float
        The `c` value that maximizes cluster separation performance.
    """

    if data.size == 0:
        raise ValueError(f"Data size is zero.")
    
    if start > end:
        raise ValueError(f"Invalid line segment [{start}, {end}).")
    
    # convert labels to integers
    lab = __process_labels(labels)

    # check if there is a label for every cell
    if len(lab) != data.shape[0]:
        # check if data has to be transposed
        # we want cells as rows and genes as columns
        raise ValueError(f"Label vector has length {len(lab)}, while the data matrix is of shape {data.shape}.")
    
    return __optimise_c_internal(data, lab, num_markers, start, end, step_size, verbose=verbose)
