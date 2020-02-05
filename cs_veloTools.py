# custom function to also reduce umap data

def filter_cells(self, bool_array: np.ndarray) -> None:
    """Filter cells using a boolean array.
    Arguments
    ---------
    bool_array: np.ndarray (size )
        array describing the cells to keep (True).
    Return
    ------
    Nothing but it removes some cells from S and U.
    """
    self.S, self.U, self.A = (X[:, bool_array] for X in (self.S, self.U, self.A))
    self.initial_cell_size = self.initial_cell_size[bool_array]
    self.initial_Ucell_size = self.initial_Ucell_size[bool_array]
    self.ca = {k: v[bool_array] for k, v in self.ca.items()}
    
    # for object that have been normalized already
    try:
        self.S_sz, self.S_norm, self.U_sz, self.U_norm = (X[:, bool_array] for X in (self.S_sz, self.S_norm, self.U_sz, self.U_norm))
        self.cell_size = self.cell_size[bool_array]
        self.Ucell_size = self.Ucell_size[bool_array]
        self.norm_factor = self.norm_factor[bool_array]    
        self.Unorm_factor = self.Unorm_factor[bool_array]    
    except:
        pass
    
    # for objects that have been dimension reduced
    try:
        self.pcs = self.pcs[bool_array]  # type: np.ndarray
    except:
        pass
    try:
        self.ts = self.ts[bool_array]  # type: np.ndarray
    except:
        pass
    try:
        self.um = self.um[bool_array]  # type: np.ndarray
    except:
        pass
    try:
        self.embedding = self.embedding[bool_array]  # type: np.ndarray
    except:
        pass
    try:
        self.delta_embedding_random = self.delta_embedding_random[bool_array]  # type: np.ndarray
    except:
        pass
    try:
        self.delta_embedding = self.delta_embedding[bool_array]  # type: np.ndarray
    except:
        pass
    
    
    
    # for objects which have been score_cv_vs_mean'd
    try:
        self.size_factor = self.size_factor[bool_array]  # type: np.ndarray
    except:
        pass
    
    # for objects that have cluster assignments
    try:
        self.cluster_labels = self.cluster_labels[bool_array]  # type: np.ndarray
        self.colorandum = self.colorandum[bool_array, :]  # type: np.ndarray
    except AttributeError:
        pass