import numpy as np
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp
import pandas

class MarkerGeneLPSolver:
    def __init_vars(self):
        # general solver parameters
        self.__data = None
        self.__labels = None
        self.__num_markers = None
        self.__s_squared = None
        self.__method = None
        self.__weighted = False

        self.__ilp = False

        # label to indices of points dictionary
        self.__idx = {}

        # label to data
        self.__cells_by_label = {}

        # total number of groups of constraints
        self.__m = None

        # list which contains all labels (once each)
        self.__unique_labels = None

        # dictionaries used to cache variance and means
        self.__variance_dict = {}
        self.__mean_dict = {}

    def __init__(self, data, labels, num_markers, s, weighted=False, k = 0.5, method="alternate", ilp=False):
        self.__init_vars()

        self.__data = data
        self.__num_markers = num_markers
        self.__s_squared = s * s
        self.__method = method
        self.__weighted = weighted
        self.__num_cells, self.__num_genes = data.shape
        self.__k = k
        self.__ilp = ilp
        self.__labels = labels

        # compute label to indices dict
        for i in range(len(self.__labels)):
            label = self.__labels[i]
            if label not in self.__idx:
                self.__idx[label] = []
            self.__idx[label].append(i)

        # compute label to data dict
        for (lab, indices) in self.__idx.items():
            self.__cells_by_label[lab] = np.asarray(self.__data[indices, :])

        # get unique labels
        self.__unique_labels = list(self.__idx.keys())

        # compute number of constraints
        self.__m = self.__get_m()

    def __get_m(self):
        # number of constraint is equal to the number of pairs
        return (len(self.__unique_labels) * (len(self.__unique_labels) - 1)) // 2
    
    def __get_data_by_label(self, label):
        return self.__cells_by_label[label]

    def __get_variance(self, label, mean):
        # first we try to get cached variance if it exists
        value = self.__variance_dict.get(label, None)
        if value is None:
            # get data for specific cell type
            data = self.__get_data_by_label(label)

            # cached variance not found - we calculate it
            n = data.shape[0]
            
            diff = np.square(data - mean)
            var = np.sum(diff, axis=0) # sum by columns

            # convert to ndarray if needed
            if isinstance(var, np.matrix):
                var = var.A.reshape(-1) 

            value = var / (n - 1)

            # cache the recently calculated variance
            self.__variance_dict[label] = value

        return value
    
    def __get_mean(self, label):
        # first we try to get cached mean if it exists
        value = self.__mean_dict.get(label, None)
        if value is None:
            # nothing is cached, we have to calculted this and cache it

            # split the dataset by labels i and j
            data = self.__get_data_by_label(label)

            value = np.mean(data, axis=0)
            self.__mean_dict[label] = value
        
        return value

    def __singleton(self, val):
        return np.array([val])
        
    def __get_weights(self):
        if self.__weighted:
            # in the weighted case, we have to calculate the number of 
            # pairwise distances between two sets of cell types
            coeffs = np.empty(self.__m)
            cnt = 0
            for i in range(len(self.__unique_labels)):
                for j in range(len(self.__unique_labels)):
                    # allow only pairs (i, j) such that i < j
                    if j <= i:
                        continue

                    c1 = len(self.__idx[self.__unique_labels[i]])
                    c2 = len(self.__idx[self.__unique_labels[j]])

                    coeffs[cnt] = np.log(c1 + c2)
                    cnt += 1
        else:
            # in the unweighted case, all weights are 1
            coeffs = np.ones(self.__m)
            
        return np.concatenate((
                   np.zeros(self.__num_genes), # 0 attached to each alpha in the objective func
                   np.zeros(self.__m), # 0 attached to each z_ij in the objective
                   coeffs
               ))

    def __get_constraints(self):
        n, d = self.__data.shape        

        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()

            with gp.Model("SepSolve", env=env) as M:
                # variables:
                # d alphas
                # m z's 
                # m betas

                # variable upper bounds:
                u = np.concatenate((
                    np.array([1 for i in range(d)]), # upper bound for alphas is 1 (or fixed value)
                    float('inf') * np.ones(2 * self.__m) # upper bound for betas and z's is inf
                ))

                # variable lower bounds:
                l = np.concatenate((
                    np.array([0 for i in range(d)]), # lower bound for alphas is 0 (or fixed value)
                    float(0.0) * np.ones(2 * self.__m) # lower bound for z's and betas is 0
                ))
                
                if self.__ilp:
                    vtypes = [GRB.BINARY] * d + [GRB.CONTINUOUS] * (2 * self.__m) # alphas are integral, z's and betas continuous
                    x = M.addMVar(shape = d + 2 * self.__m, vtype = vtypes)
                else:
                    x = M.addMVar(shape = d + 2 * self.__m, vtype = GRB.CONTINUOUS, ub = u, lb = l)

                # coefficients in the objective function - minimise the sum of betas
                c = self.__get_weights()
                M.setObjective(c @ x, GRB.MINIMIZE)

                # generate (or load a cached) system matrix
                A, b = self.__get_system_matrix()

                # add constraints
                M.addConstr(A @ x >= b)

                # optimize
                M.optimize()             

                obj = M.ObjVal
                x = (M.x)[:d] # get alphas
                y = (M.x)[-self.__m:] # get betas

        return x, y, obj
    
    def __get_system_matrix(self):
        n, d = self.__data.shape

        # fill matrix A with numbers
        all_vals, all_rows, all_cols, rhs = self.__get_cons_matrix()

        # finally, add the dimension constraints
        # -sum alpha >= -num_markers
        all_vals.append(-np.ones(d)) # d ones attached to alphas
        all_cols.append(np.arange(d)) # d columns mathing d alphas
        all_rows.append(np.array([3 * self.__m] * d)) # all d values are in the final row

        # sum alpha >= num_markers
        all_vals.append(np.ones(d)) # d ones attached to alphas
        all_cols.append(np.arange(d)) # d columns mathing d alphas
        all_rows.append(np.array([3 * self.__m + 1] * d)) # all d values are in the final row

        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        vals = np.concatenate(all_vals)
        
        # we're done! let's generate the matrix
        A = sp.csr_matrix(
            (vals, (rows, cols)),
            shape = (3 * self.__m + 2, d + 2 * self.__m) # total 3m constraints and d + 2m vars, 2 extra cons for maximum dimension of alpha
        )
        
        # right hand side
        b = np.concatenate((
            rhs,
            np.array([-int(self.__num_markers)]), # final dimensionality constraint
            np.array([int(self.__num_markers)]), # final dimensionality constraint
        ))

        return A, b
    
    def __get_cons_matrix(self):           
        return self.__get_cons_matrix_alternate()
    
    def __get_cons_matrix_alternate(self):
        n, d = self.__data.shape
        
        all_vals = []
        all_rows = []
        all_cols = []
        all_b = []

        cnt = 0
        # iterate through all labels
        # NOTE: we assume that the labels are all integers from 0 to k
        for i in range(len(self.__unique_labels)):
            for j in range(len(self.__unique_labels)):
                # allow only pairs (i, j) such that i < j
                if j <= i:
                    continue 

                # calculating menas can be slow 
                # this function will cache and re-use them
                m1 = self.__get_mean(i)
                m2 = self.__get_mean(j)
                diff = np.square(m1 - m2)

                # convert to ndarray if needed - operations on sparse matrices return np.matrix
                # this messes up concatenation later on
                if isinstance(diff, np.matrix):
                    diff = diff.A.reshape(-1)

                # calculating variance is expensive
                # this function will cache and re-use them
                s1 = self.__get_variance(i, m1)
                s2 = self.__get_variance(j, m2)

                s1_sum = np.sum(s1)
                s2_sum = np.sum(s2)

                # for this pair we generate 3 constraints:
                values = np.concatenate((
                    # first constraint
                    diff, # first d alphas
                    self.__singleton(float(-1.0)), # coeff with z_cnt

                    # second constraint
                    -(self.__s_squared / 2) * s1, # d values paired with alphas
                    self.__singleton(float(1.0)), # coef with z_cnt
                    self.__singleton(s1_sum * self.__k), # coef with beta_cnt

                    # third constraint
                    -(self.__s_squared / 2) * s2, # d values paired with alphas
                    self.__singleton(float(1.0)), # coef with z_cnt
                    self.__singleton(s2_sum * self.__k), # coef with beta_cnt
                ))

                rows = np.concatenate((
                    # first constraint
                    np.array([3 * cnt] * (d + 1)), # d + 1 values in first row: d alphas and one z

                    # second constraint
                    np.array([3 * cnt + 1] * (d + 2)), # d + 1 values in second row: d alphas, z and beta

                    # third constraint
                    np.array([3 * cnt + 2] * (d + 2)), # d + 1 values in third row: d alphas, z and beta
                ))

                cols = np.concatenate((
                    # first constraint
                    np.arange(d), # first d columns are with alpha
                    self.__singleton(d + cnt), # column for z_cnt variable

                    # second constraint
                    np.arange(d), # first d columns are with alpha
                    self.__singleton(d + cnt), # column for z_cnt variable
                    self.__singleton(d + self.__m + cnt), # column for beta_cnt variable

                    # third constraint
                    np.arange(d), # first d columns are with alpha
                    self.__singleton(d + cnt), # column for z_cnt variable
                    self.__singleton(d + self.__m + cnt), # column for beta_cnt variable
                ))

                # right hand side for these three constraints
                b = np.array([
                    float(0.0), # first constraint  
                    self.__k * self.__s_squared * s1_sum / 2, # second constraint
                    self.__k * self.__s_squared * s2_sum / 2, # third constraint
                ])

                all_vals.append(values)
                all_rows.append(rows)
                all_cols.append(cols)
                all_b.append(b)
                cnt += 1

        return all_vals, all_rows, all_cols, np.concatenate(all_b)
    
    def ranking(self, x):
        # performs simple ranking on the solution
        return sorted(range(len(x)), key=lambda i: x[i], reverse=True)[: self.__num_markers]

    def Solve(self):
        x, y, obj = self.__get_constraints()
        return x, y, obj
  
def __get_markers__internal(data, labels, num_markers, s=0.4, k=0.5, weighted=False, ilp=False):
    solver = MarkerGeneLPSolver(data, labels, num_markers, s, weighted=weighted, k=k, ilp=ilp)
    x, _, _ = solver.Solve()
    x = solver.ranking(x)
    return x

def __process_labels(labels):
    # update labels
    if isinstance(labels, pandas.core.series.Series):
        return labels.astype('category').astype("category").cat.codes.to_numpy()
        
def get_markers(data, labels, num_markers, c=0.4, ilp=False):
    lab = __process_labels(labels)
    return __get_markers__internal(data, lab, num_markers, s=c, ilp=ilp)