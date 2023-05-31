import numpy as np
import sys


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
    The first layer is guaranteed to have at most {n_fl_co} nodes.
    This method ensures each layer is tightly packed on its diameter, except the outermost one.
    
    This has the potential to be expanded on where n_fl_co is calculated with the outermost layer's
    spacing in mind. If that layer's spacing is too sparse, then n_fl_co can be iterated up or down.
    
    Returns:
        - R_arr: array of length {layers} for each radius
        - n_arr: array of length {layers} for number of nodes per layer
    '''
    n_fl = min(n_tot, n_fl_co)

    R_arr = [big_r(n_fl, r)]

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


def get_xy(n, n_fl_co=20, r=1050):
    n_tot = int(n)
    n_fl_co = int(n_fl_co)
    r = float(r)

    assert n > 0, f"{n}: Total number of nodes must be an integer > 0"
    assert n_fl_co > 0, f"{n_fl_co}: At least 1 node in the first layer is required"
    assert r > 0, f"{r}: r must be a float > 0"

    R_arr, n_arr = even_spacing_method(n_tot, n_fl_co, r)
    Xs, Ys = node_coords(R_arr, n_arr)

    return Xs, Ys, R_arr, n_arr
