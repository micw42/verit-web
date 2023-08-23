import numpy as np
import networkx as nx
from numpy.lib.stride_tricks import sliding_window_view
import math
import pandas as pd


def get_ring_coord(n, R, offset=False):
    '''
    Place {n} nodes evenly along a circle of radius {R}.
    Arguments:
        - n (int): number of nodes
        - R (float): radius
        - offset (bool): shift angle by half a unit
    Returns:
        - x, y: two vectors -- one for the x position and one for the y
    '''
    theta = np.linspace(0, 2*np.pi - (2*np.pi/n), n)
    if offset:
        theta = theta + (np.pi/n)

    x = R * (np.sin(theta))
    y = R * (np.cos(theta))

    return x, y


# Function to calculate the ring's radius given padding, number of nodes, and radius of each node
big_r = lambda n, r, p=0: n * (2*r + p) / (2*np.pi)

def even_number_method(n, r, layers=1):
    '''
    Get radius values for each layer given {n} nodes with {r} radius each.
    This method ensures the same number of nodes are in each layer.
    
    Returns:
        - R_arr: array of length {layers} for each radius
        - n_arr: array of length {layers} for number of nodes per layer
    '''
    n_arr = np.array([len(x) for x in np.array_split(np.zeros(n), layers)])
    n_arr = np.flip(n_arr)    # Flip so the more remainders are outside
    
    # Initialize first radius where padding is 0
    # Add 2r/(n-1) to handle cases where {n} is small
    R_arr = np.repeat(big_r(n_arr[0], r) + 2*r/(n-1), layers)
    # The next radii go outwards linearly
    R_arr = np.array([x + i*(2*r) for i, x in enumerate(R_arr)])
    
    return R_arr, n_arr


def even_spacing_method(n_tot, n_fl_co, r):
    '''
    Get radius values for each layer given {n_tot} nodes with {r} radius each. 
    The first layer is guaranteed to have at most {n_fl_co} nodes (number of first-layer cutoff).
    This method ensures each layer is tightly packed on its diameter, except the outermost one.
    
    This has the potential to be expanded on where n_fl_co is calculated with the outermost layer's
    spacing in mind. If that layer's spacing is too sparse, then n_fl_co can be iterated up or down.
    
    Returns:
        - R_arr: array of length {layers} for each radius
        - n_arr: array of length {layers} for number of nodes per layer
    '''
    n_fl = min(n_tot, n_fl_co)

    R_arr = [max(big_r(n_fl, r), 4*r)]

    n_arr = [n_fl]
    n_fitted = n_fl
    while n_fitted < n_tot:
        R_i = R_arr[-1] + 2*r

        outer_diam = R_i * 2 * np.pi
        i_fitted = outer_diam // (2*r)

        n_fitted += i_fitted
        n_arr.append(i_fitted)
        R_arr.append(R_i)

    # Last one needs special treatment
    n_arr[-1] = n_arr[-1] - (np.sum(n_arr) - n_tot)
    n_arr = [int(x) for x in n_arr]

    return R_arr, n_arr


def node_coords(R_arr, n_arr):
    '''
    Assign coordinate values for an array of {R_arr} ring radii -- each ring containing {n_arr[i]} nodes.
    Returns:
        - Xs, Ys: arrays of X and Y values each with length {n}
    '''

    Xs = []
    Ys = []
    offset = False
    for r, ns in zip(R_arr, n_arr):
        X, Y = get_ring_coord(ns, r, offset=offset)
        Xs.append(X)
        Ys.append(Y)
        offset = not offset

    Xs = np.concatenate(Xs)
    Ys = np.concatenate(Ys)
    
    return Xs, Ys


def get_xy(n, n_fl_co=20, r=1050, offset_x=0, offset_y=0):
    n_tot = int(n)
    n_fl_co = int(n_fl_co)
    r = float(r)

    assert n >= 0, f"{n}: Total number of nodes must be an integer > 0"
    assert n_fl_co > 0, f"{n_fl_co}: At least 1 node in the first layer is required"
    assert r > 0, f"{r}: r must be a float > 0"
    
    if n == 0:
        return None, None, None, None

    R_arr, n_arr = even_spacing_method(n_tot, n_fl_co, r)
    Xs, Ys = node_coords(R_arr, n_arr)
    
    Xs = Xs + offset_x
    Ys = Ys + offset_y

    return Xs, Ys, R_arr, n_arr

def SQ_layered_concentric(qnodes_df, qedges_df):
     # Sort the nodes by thickness to be arranged polarly
    query_id = qnodes_df[qnodes_df["depth"] == 0].iloc[0]["Id"]

    nq1 = qedges_df[qedges_df.source == query_id][["target", "thickness"]].rename(columns={"target": "Id"})
    nq2 = qedges_df[qedges_df.target == query_id][["source", "thickness"]].rename(columns={"source": "Id"})

    nq_df = pd.concat([nq1, nq2])
    nq_df = nq_df.groupby("Id").max().reset_index()
    qnodes_df = qnodes_df.merge(nq_df, on="Id", how="left")
    qnodes_df = qnodes_df.sort_values(["depth", "thickness"], ascending=[True, False])

    # Calculate x, y coordinates for each node
    Xs, Ys, _, __ = get_xy(len(qnodes_df)-1, n_fl_co=20, r=1050)
    Xs = [0] + list(Xs)    # First index is 0 because it's the query node
    Ys = [0] + list(Ys)
    
    qnodes_df["lc_X"] = Xs; qnodes_df["lc_Y"] = Ys
    
    return qnodes_df
    

def layered_concentric(qnodes_df):
    qnodes_df = qnodes_df.sort_values("Type", ascending=False)
    vc_nodes = qnodes_df.Type.value_counts()

    n_query = vc_nodes["Query"]
    n_links = vc_nodes.sum() - vc_nodes["Query"]

    # Compute X and Y for concentric layout
    Xs = []; Ys = []

    r1 = 500
    Xs_q, Ys_q, R_arr_q, n_arr_q = get_xy(n_query, r=r1)
    Xs.extend(Xs_q); Ys.extend(Ys_q)

    if n_links>0:
        r2 = 100
        n_fl_co_d = 2 * np.pi * (R_arr_q[-1] + 3*r1) / (2 * r2)
        Xs_d, Ys_d, R_arr_d, n_arr_d = get_xy(n_links, n_fl_co_d, r=r2)
        Xs.extend(Xs_d); Ys.extend(Ys_d)
    
    qnodes_df["lc_X"] = Xs; qnodes_df["lc_Y"] = Ys
    
    return qnodes_df


def assign_cluster(qedges_df, qnodes_df):
    ''' Assign each node to a query "cluster"
    Argument:
        qedges_df (pd.DataFrame): edge dataframe for query
        qnodes_df (pd.DataFrame): node dataframe for query
    Returns:
        qnodes_df (pd.DataFrame): updated nodes dataframe with extra "clust" column
        clust_sizes(pd.Series): size of each query cluster
    '''
    qnodes_df = qnodes_df.reset_index(drop=True)
    qedges_df = qedges_df.reset_index(drop=True)
    
    # Prepare the cluster column. Note: query nodes will have their own Id as their cluster.
    qnodes_df["clust"] = qnodes_df["Id"]

    qG = nx.from_pandas_edgelist(qedges_df, edge_attr=True, source="source", target="target", create_using=nx.Graph())
    adj_mat = nx.adjacency_matrix(qG, nodelist=qnodes_df["Id"], weight="thickness")

    # Get query and non-query nodes, construct query-non-query adjacency matrix
    queries = qnodes_df[qnodes_df["Type"] == "Query"]["Id"]
    q_i = np.array(queries.index)
    nq_i = np.setdiff1d(range(len(qnodes_df)), q_i)
    adj_mat = adj_mat[np.ix_(nq_i, q_i)]

    # Divide columns by column sums to normalize over each query node's baseline "popularity"
    colsum = np.sum(adj_mat, axis=0)
    adj_mat = np.nan_to_num(adj_mat/colsum, 0)

    # Assign each non-query node to a query node
    assign = np.asarray(np.argmax(adj_mat, axis=1).squeeze()).squeeze()
    assign = queries.reset_index(drop=True)[assign]    # Now want to index specific to the list of queries
    assign.index = nq_i
    qnodes_df.loc[nq_i,"clust"] = assign
    
    clust_sizes = qnodes_df.groupby("clust").size().sort_values(ascending=False)
    
    # Sort nodes table to be in the order of assigning X and Y coordinates
    qnodes_df.clust = qnodes_df.clust.astype("category")
    qnodes_df.clust = qnodes_df.clust.cat.set_categories(clust_sizes.index)
    qnodes_df = qnodes_df.sort_values(["clust", "Type"], ascending=[True, False]).reset_index(drop=True)
    
    return qnodes_df, clust_sizes


def cluster_xy(ring_R, Rs, offset=False, arr_family=False, denom=4):
    ''' Compute the center location of clusters given the radius of the current layer and variable cluster radii
    Arguments:
        ring_R (float): radius of the layer that the cluster should reside along
        Rs (arr): array of cluster radii
        offset (bool): stagger X and Y values to help visualization
        arr_family (bool): arrange clusters by family (avoids monotonic decreasing cluster size)
        denom (int): consider 1/{denom} of the total clusters as big (only applicable if gini coefficient not above threshold)
    Returns:
        cX, cY (arr): coordinates for each cluster at current layer
    '''
    # Calculate the circumference
    circ = ring_R*2*np.pi
    
    # Calculate number of nodes counted as big in family assignment
    gc = gini_coefficient(Rs)
    if gc > 0.3:
        n_big = np.sum(Rs > np.mean(Rs))
    else:
        n_big = math.floor(max(denom, len(Rs))/denom) #minimum of 1 big node
    
    #Get index for family-sorted Rs
    if arr_family:
        family = assign_family(Rs, n_big=n_big)
        family_idx = np.argsort(family)
    else: 
        family_idx = np.arange(len(Rs))

    # Sort Rs by family
    Rs = Rs[family_idx]
    
    # The cumulative sum of sliding window sum (of radii) corresponds to the center locations along a linearized ring
    # Handle if there is only 1 R in Rs
    if len(Rs) == 1:
        Rs = np.insert(Rs, -1, 0)
    lin_centers = np.insert(np.cumsum(np.sum(sliding_window_view(Rs, window_shape = 2), axis = 1)), 0, 0)
    
    # Compute the remainder space within the layer and evenly distribute it as padding between clusters
    remainder = circ - (lin_centers[-1]+Rs[-1]+Rs[0])
    assert remainder >= 0, f"Layer overcrowded: remainder of {remainder}"   # Check if remainder is negative
    remainder_pad = remainder / (len(Rs))
    
    # Push the centers by the remainder padding. Automatically accounts for inter-cluster padding.
    lin_centers = [lin_centers[i] + remainder_pad*i for i in range(len(lin_centers))]
    lin_centers = np.array(lin_centers) / circ
    lin_centers = lin_centers*2*np.pi
    
    # Restore original order
    lin_centers = lin_centers[np.argsort(-Rs)]
    
    # Stagger adjacent layers for readability (this is not the optimal way to do it)
    if not offset:
        cX = ring_R * np.cos(lin_centers)
        cY = ring_R * np.sin(lin_centers)
    else:
        cX = ring_R * np.cos(lin_centers+(np.pi/2))
        cY = -ring_R * np.sin(lin_centers+(np.pi/2))
    
    return cX, cY


def cluster_layer(clust_sizes, icp=10000, n_fl_co=20, r=1050, arr_family=False, denom=4):
    ''' Compute the coordinates for each cluster
    Arguments:
        clust_sizes (arr): how many nodes in each cluster
        icp (float): inter-cluster padding. How much each cluster should be spaced out
        **kwargs: for even_spacing_method (so, radius and number of first-layer cutoff)
    
    Returns:
        cXs, cYs (arr): coordinates for each cluster
    '''
    assert (np.array(clust_sizes) == -np.sort(-clust_sizes)).all(), "Sort by descending first"
      
    # First retrieve each cluster's radius
    R_maxes = []
    for n in clust_sizes:
        R_arr, n_arr = even_spacing_method(n, n_fl_co=n_fl_co, r=r)
        R_maxes.append(max(R_arr))
    R_maxes = np.array(R_maxes)

    # Initialize cluster membership and radii
    layer_assign = np.zeros(len(R_maxes))
    layer_assign[0] = 1
    R_arr = [0, R_maxes[0]+icp+R_maxes[1]]

    # Cumulative sum is used to track how many variable-radius nodes can fit in each layer
    layer = 2
    csum = np.cumsum((2*R_maxes) + icp)
    csum = csum - csum[0]
    
    cXs = [0]
    cYs = [0]

    offset = False    # Used to improve visibility of adjacent layers
    # While at least one cluster has not been assigned yet...
    while np.product(layer_assign) == 0:
        curr_circ = 2*np.pi*R_arr[-1]

        # Greater than 0 because csum will be subtracted by the csum immediately exceeding current layer's capacity
        i_fit = (csum <= curr_circ) & (csum > 0)
        # Assign layer to the clusters that fit
        layer_assign[i_fit] = layer
        
        # Record the maximum radius fitting in the current layer to update the next layer's radius
        R_max_layer = max(R_maxes[np.where(i_fit)])
        
        # Get adjustment coordinate of cluster
        cX, cY = cluster_xy(R_arr[-1], R_maxes[i_fit], offset=offset, arr_family=arr_family, denom=denom)
        offset = not offset
        cXs.extend(cX); cYs.extend(cY)
        
        R_arr.append(R_arr[-1]+icp+R_max_layer)
        
        # Dummy way of avoiding error when loop probably met exit condition
        try:
            i_fit_max = max(np.squeeze(np.where(i_fit)))
            csum = csum - csum[i_fit_max]    # Update csum to account for fitting the previous layer
        except:
            break
        
        layer += 1
    
    return cXs, cYs

def gini_coefficient(x):
    """Compute Gini coefficient of array of values"""
    diffsum = 0
    for i, xi in enumerate(x[:-1], 1):
        diffsum += np.sum(np.abs(xi - x[i:]))
    return diffsum / (len(x)**2 * np.mean(x))

def assign_family(clust_sizes, n_big=3):
    ''' Assigns each cluster to a family
    Arguments:
        clust_sizes (arr): how many nodes in each cluster, sorted increasing to decreasing
        n_big (int): top n nodes to consider as "big"
    
    Returns:
        family (arr): int representing the family each cluster belongs to
    '''
    n_big = min(n_big, len(clust_sizes))
    n_tile = math.ceil(len(clust_sizes)/n_big)
    family = np.tile(np.arange(n_big), n_tile)[:len(clust_sizes)]
    return family

def cluster_layered_concentric(qnodes_df, qedges_df, r=100, icp=500, arr_family=False, denom=4):
    qnodes_df, clust_sizes = assign_cluster(qedges_df, qnodes_df)
    cXs, cYs = cluster_layer(clust_sizes-1, r=r, icp=icp, arr_family=arr_family, denom=denom)

    full_Xs = []; full_Ys = []
    for i in range(len(clust_sizes)):
        # First node is query node and center of cluster
        Xs, Ys, _, _ = get_xy(clust_sizes[i]-1, offset_x=cXs[i], offset_y=cYs[i], r=r)

        full_Xs.append(cXs[i]); full_Ys.append(cYs[i])
        if Xs is not None:
            full_Xs.extend(Xs); full_Ys.extend(Ys)
    
    qnodes_df["clc_X"] = full_Xs; qnodes_df["clc_Y"] = full_Ys
    
    return qnodes_df