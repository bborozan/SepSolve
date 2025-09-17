import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score

from .sepsolve_base import MarkerGeneLPSolver
from .sepsolve_fixed import SepSolveFixed
        
def __optimise_c_internal(data, labels, num_markers, start, end, step_size, verbose=False):
    # split original dataset
    data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.3) # 70/30 split

    base = MarkerGeneLPSolver(data_train, labels_train, num_markers, ilp=False)

    best_c = None
    best_score = -1

    for c in np.arange(start, end + 1e-9, step_size):
        if verbose:
            print(f"Testing separation constant {c}")

        # run SepSolve with c and get markers
        solver = SepSolveFixed(base, c)

        x, betas, obj = base.Solve(solver)
        markers = base.ranking(x)

        # copy data to avoid using slices
        data_train_markers = data_train[:, markers].copy()
        data_test_markers = data_test[:, markers].copy()

        # train logistic regression on training set
        clf = LogisticRegression(max_iter=1200, verbose=0, n_jobs=-1).fit(data_train_markers, labels_train)
        pred = clf.predict(data_test_markers)
        
        # compute macro average F1 score
        f1 = f1_score(labels_test, pred, average="macro") 

        # check if its better
        if f1 > best_score:
            best_score = f1
            best_c = c
    
    return best_c
