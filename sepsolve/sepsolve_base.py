import numpy as np
from scipy.sparse import issparse

class MarkerGeneLPSolver(): 
    def __init_vars(self):
        # general solver parameters
        self.__data = None
        self.__labels = None
        self.__num_markers = None

        self.__ilp = False

        # label to indices of points dictionary
        self.__idx = {}

        # label to data
        self.__cells_by_label = {}

        # list which contains all labels (once each)
        self.__unique_labels = None
        self.__num_clusters = 0

        # dictionaries used to cache variance and means
        self.__variance_dict = {}
        self.__mean_dict = {}


    def __init__(self, data, labels, num_markers, ilp=False):
        self.__init_vars()

        self.__data = data
        self.__num_markers = num_markers
        self.__num_cells, self.__num_genes = data.shape
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
            label_data = self.__data[indices, :]

            if issparse(label_data):
                # convert from sparse to dense
                # np.mean does not work well with sparse format
                label_data = label_data.toarray()
            else:
                # it already was a numpy array
                label_data = np.asarray(label_data)

            if label_data.size > 0:
                self.__cells_by_label[lab] = label_data

        # get unique labels
        self.__unique_labels = list(self.__idx.keys())
        self.__num_clusters = len(self.__unique_labels)

    @property
    def data(self):
        return self.__data
    
    @property
    def labels(self):
        return self.__labels

    @property
    def num_markers(self):
        return self.__num_markers
    
    @property
    def num_genes(self):
        return self.__num_genes
    
    @property
    def num_clusters(self):
        return self.__num_clusters
    
    @property
    def num_cells(self):
        return self.__num_cells
    
    @property
    def ilp(self):
        return self.__ilp
    
    @property
    def unique_labels(self):
        return self.unique_labels

    def __get_data_by_label(self, label):
        return self.__cells_by_label[label]

    def get_variance(self, label, mean):
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

            if n == 1:
                value = np.zeros(var.shape)
            else:
                value = var / (n - 1)

            # cache the recently calculated variance
            self.__variance_dict[label] = value

        return value
    
    def get_mean(self, label):
        # first we try to get cached mean if it exists
        value = self.__mean_dict.get(label, None)
        if value is None:
            # nothing is cached, we have to calculted this and cache it

            # split the dataset by labels i and j
            data = self.__get_data_by_label(label)

            value = np.mean(data, axis=0)
            self.__mean_dict[label] = value
        
        return value
    
    def ranking(self, x):
        # performs simple ranking on the solution
        return sorted(range(len(x)), key=lambda i: x[i], reverse=True)[: self.__num_markers]
    
    def Solve(self, solver):
        return solver.Solve()
    
class MarkerGeneSolver():
    def singleton(self, val):
        return np.array([val])

    def Solve():
        raise NotImplementedError()
