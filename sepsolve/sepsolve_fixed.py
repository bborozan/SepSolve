import numpy as np
import gurobipy as gp
from gurobipy import GRB

from .sepsolve_base import MarkerGeneLPSolver
from .sepsolve_base import MarkerGeneSolver 

class SepSolveFixed(MarkerGeneSolver):
    # SepSolve with a fixed choice of parameters c 
    def __init__(self, solver : MarkerGeneLPSolver, s):
        self.parent = solver
        self.__m = self.get_m()
        self.__s_squared = s * s

    def get_m(self):
        # number of constraint is equal to the number of pairs
        return (self.parent.num_clusters * (self.parent.num_clusters - 1)) // 2

    def get_weights(self):
        return np.concatenate((
            np.zeros(self.parent.num_genes), # alphas
            np.ones(self.__m) # betas
        ))
    
    def __setup_variables(self, M):
        d = self.parent.num_genes
        # variable upper bounds:
        u = np.concatenate((
            np.ones(d), # upper bound for alphas is 1 (or fixed value)
            np.inf * np.ones(self.__m) # upper bound for betas
        ))

        # variable lower bounds:
        l = np.concatenate((
            np.zeros(d), # lower bound for alphas is 0
            np.zeros(self.__m) # lower bound for betas is 0
        ))
        
        if self.parent.ilp:
            vtypes = [GRB.BINARY] * d + [GRB.CONTINUOUS] * self.__m # alphas are integral, betas continuous
            x = M.addMVar(shape = d + self.__m, vtype = vtypes)
        else:
            x = M.addMVar(shape = d + self.__m, vtype = GRB.CONTINUOUS, ub = u, lb = l)

        return x

    def __add_contraint_single(self, M, x, lhs, rhs):
        M.addConstr(lhs @ x >= rhs)

    def __add_dim_constraint(self, M, x):
        lhs = np.concatenate((np.ones(self.parent.num_genes), np.zeros(self.__m)))
        M.addConstr(lhs @ x == self.parent.num_markers)

    def __add_cons(self, M, x):
        cnt = 0
        for i in range(self.parent.num_clusters):
            m1 = self.parent.get_mean(i)
            s1 = self.parent.get_variance(i, m1)
            s1_sum = np.sum(s1)

            for j in range(self.parent.num_clusters):
                # allow only pairs (i, j) such that i < j
                if j <= i:
                    continue 

                m2 = self.parent.get_mean(j)
                diff = np.square(m1 - m2)

                # convert to ndarray if needed - operations on sparse matrices return np.matrix
                # this messes up concatenation later on
                if isinstance(diff, np.matrix):
                    diff = diff.A.reshape(-1)

                s2 = self.parent.get_variance(j, m2)
                s2_sum = np.sum(s2)

                beta1 = np.zeros(self.__m)
                beta1[cnt] = s1_sum * 0.5

                beta2 = np.zeros(self.__m)
                beta2[cnt] = s2_sum * 0.5

                lhs1 = np.concatenate((
                    diff - (self.__s_squared * 0.5) * s1,
                    beta1
                ))
                rhs1 = (self.__s_squared * 0.25) * s1_sum

                lhs2 = np.concatenate((
                    diff - (self.__s_squared * 0.5) * s2,
                    beta2
                ))
                rhs2 = (self.__s_squared * 0.25) * s2_sum

                self.__add_contraint_single(M, x, lhs1, rhs1)
                self.__add_contraint_single(M, x, lhs2, rhs2)

                cnt += 1

            self.__add_dim_constraint(M, x)
        
    def Solve(self):
        d = self.parent.num_genes     

        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()

            if self.__s_squared is not None and self.__s_squared > 0:
                with gp.Model("SepSolve", env=env) as M:
                    # variables:
                    # d alphas
                    # m betas
                    x = self.__setup_variables(M)

                    # coefficients in the objective function - minimise the sum of betas
                    c = self.get_weights()
                    M.setObjective(c @ x, GRB.MINIMIZE)

                    self.__add_cons(M, x)

                    # optimize
                    M.optimize()              
                    obj = M.ObjVal
                    x = (M.x)[:d] # get alphas
                    y = (M.x)[-self.__m:] # get betas
                    return x, y, obj
                