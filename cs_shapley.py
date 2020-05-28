def shapley_values(G, k=None):
    # modified from betweenness centrality from network x
    # TODO : fix description
    # TODO : add second term
    r"""Compute the shapley centrality for nodes.

    Betweenness centrality of a node $v$ is the sum of the
    fraction of all-pairs shortest paths that pass through $v$

    .. math::

       c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}

    where $V$ is the set of nodes, $\sigma(s, t)$ is the number of
    shortest $(s, t)$-paths,  and $\sigma(s, t|v)$ is the number of
    those paths  passing through some  node $v$ other than $s, t$.
    If $s = t$, $\sigma(s, t) = 1$, and if $v \in {s, t}$,
    $\sigma(s, t|v) = 0$ [2]_.

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    k : int, optional (default=None)
      If k is not None use k node samples to estimate betweenness.
      The value of k <= n where n is the number of nodes in the graph.
      Higher values give better approximation.

    normalized : bool, optional
      If True the betweenness values are normalized by `2/((n-1)(n-2))`
      for graphs, and `1/((n-1)(n-2))` for directed graphs where `n`
      is the number of nodes in G.

    weight : None or string, optional (default=None)
      If None, all edge weights are considered equal.
      Otherwise holds the name of the edge attribute used as weight.

    endpoints : bool, optional
      If True include the endpoints in the shortest path counts.

    Returns
    -------
    nodes : dictionary
       Dictionary of nodes with betweenness centrality as the value.

    See Also
    --------
    edge_betweenness_centrality
    load_centrality

    Notes
    -----
    The algorithm is from Ulrik Brandes [1]_.
    See [4]_ for the original first published version and [2]_ for details on
    algorithms for variations and related metrics.

    For approximate betweenness calculations set k=#samples to use
    k nodes ("pivots") to estimate the betweenness values. For an estimate
    of the number of pivots needed see [3]_.

    For weighted graphs the edge weights must be greater than zero.
    Zero edge weights can produce an infinite number of equal length
    paths between pairs of nodes.

    References
    ----------
    .. [1] Ulrik Brandes:
       A Faster Algorithm for Betweenness Centrality.
       Journal of Mathematical Sociology 25(2):163-177, 2001.
       http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
    .. [2] Ulrik Brandes:
       On Variants of Shortest-Path Betweenness
       Centrality and their Generic Computation.
       Social Networks 30(2):136-145, 2008.
       http://www.inf.uni-konstanz.de/algo/publications/b-vspbc-08.pdf
    .. [3] Ulrik Brandes and Christian Pich:
       Centrality Estimation in Large Networks.
       International Journal of Bifurcation and Chaos 17(7):2303-2318, 2007.
       http://www.inf.uni-konstanz.de/algo/publications/bp-celn-06.pdf
    .. [4] Linton C. Freeman:
       A set of measures of centrality based on betweenness.
       Sociometry 40: 35–41, 1977
       http://moreno.ss.uci.edu/23.pdf
    """
    shapley_values = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    nodes = G
    
    for s in nodes:
        # single source shortest paths
        S, P, sigma, D = _single_source_shortest_path_basic(G, s)
        # accumulation
        shapley_values = _accumulate_basic(shapley_values, S, P, sigma, D, s)
    # rescaling
    shapley_values = _rescale(shapley_values, len(G), normalized=True,
                           directed=G.is_directed(), k=k)
    return shapley_values

def _single_source_shortest_path_basic(G, s):
    S = []   # list of nodes, in order
    P = {}   # dictionary of nodes to list of ?
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # dictionary of target nodes and # of shortest paths to target node from s
    D = {}   # dictionary of nodes to distances
    
    sigma[s] = 1.0    # change the value of sigma of current node (s) to 1
    D[s] = 0          # let the current distance be 0
    Q = [s]           # Q is the list of origins (which we only allow 1)
    
    while Q:   # use BFS to find shortest paths
        v = Q.pop(0)  # the current start node
        S.append(v)   # add the node that we have searched from
        Dv = D[v]     # the current distance of the current node (s)
        sigmav = sigma[v]   # the current sigma of the current node (s)
        
        for w in G[v]:      # for every node in G connected to v
            if w not in D:  # if w is not yet in D, add it to Q 
                Q.append(w)
                D[w] = Dv + 1  # and the distance between v and w is at least 1 larger than what we had before
            if D[w] == Dv + 1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v)  # predecessors (where one must go to get to target?)
    return S, P, sigma, D

def _accumulate_basic(shapleyness, S, P, sigma, D, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            if D[v] != 0:
                delta[v] += sigma[v] * coeff / D[v]
        if w != s:
            shapleyness[w] += delta[w]
    return shapleyness

def _rescale(betweenness, n, normalized, directed=False, k=None):
    if normalized:
        if n <= 2:
            scale = None  # no normalization b=0 for all nodes
        else:
            scale = 1.0 / ((n - 1) * (n - 2))
    else:  # rescale by 2 for undirected graphs
        if not directed:
            scale = 0.5
        else:
            scale = None
    if scale is not None:
        if k is not None:
            scale = scale * n / k
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness