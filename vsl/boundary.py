import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from sklearn.neighbors import KDTree
def plot_circle(centers,rs,ax):
    N = centers.shape[0]
    for i in range(N):
        theta = np.arange(0, 2*np.pi, 0.01)
        x = centers[i,0] + rs[i] * np.cos(theta)
        y = centers[i,1] + rs[i] * np.sin(theta)
        ax.plot(x, y, 'b-', alpha=0.1)

def edge_check_vaild(e,tree,r,err):
    xp = e[0]
    xq = e[1]
    L = np.sqrt(np.dot(xq-xp,xq-xp))
    if L > 2*r:
        return False, -1
    vec = (xq-xp)/L# the vector from p to q
    normal = np.array([vec[1], -vec[0]])
    c1 = (xp + xq) / 2 + normal * np.sqrt(r**2-(L/2)**2)
    c2 = (xp + xq) / 2 - normal * np.sqrt(r**2-(L/2)**2)
    c = np.array([c1, c2])
    count = tree.query_radius(c, r=r+err, return_distance=False, count_only=True, sort_results=False)
    if count[0] <= 2:
        return True, c[0]
    elif count[1] <= 2:
        return True, c[1]
    else:
        return False, -1


def boundary_extract(points, alpha, err=10e-3):
    """
    Here, parameter err was place, because there are errors when calculating distance
    meanwhile, this err was different for different scaling 2D point cloud
    so, a parameter was placed here to considering the calculation errors
    """
    R = 1 / alpha
    pts = np.copy(points)
    tree = KDTree(pts, leaf_size=2)
    tri = Delaunay(pts)
    s = tri.simplices
    N = s.shape[0]
    i = 0
    edges = []
    centers = []
    while i <= N - 1:
        if s[i, 0] == -1:
            i = i + 1
            continue
        p3 = s[i]
        e1 = np.array([points[p3[0], :], points[p3[1], :]])
        e2 = np.array([points[p3[1], :], points[p3[2], :]])
        e3 = np.array([points[p3[0], :], points[p3[2], :]])
        e = [e1, e2, e3]
        for j in range(3):
            flag, center = edge_check_vaild(e[j], tree, R, err)
            if flag:
                edges.append(e[j])
                centers.append(center)
        nb = tri.neighbors[i]
        nb_valid = nb[nb != -1]
        #nb_valid_num = nb_valid.shape[0]
        #s[nb_valid] = -1
        i = i + 1
    return edges, centers


def show_edge(edges, ax, label='SCSP', z=None, color='red', linewidth=1, alpha=0.8):
    import networkx as nx
    edf = np.array(edges)
    visit = pd.DataFrame(np.zeros(len(edf)), columns=['visited'], dtype=bool)
    edf1 = pd.DataFrame(edf[:,0], columns=['x1', 'y1'])
    edf2 = pd.DataFrame(edf[:, 1], columns=['x2', 'y2'])
    edfs = pd.concat([edf1, edf2, visit], axis=1)
    edfs = edfs.drop_duplicates()
    G = nx.Graph()
    edges0 = [((edfs.iloc[i, 0], edfs.iloc[i, 1]),
               (edfs.iloc[i, 2], edfs.iloc[i, 3])) for i in range(len(edfs))]
    G.add_edges_from(edges0)

    # draw the circle
    lines = nx.cycle_basis(G)
    for line in lines:
        x = [p[0] for p in line]# add the first vertex to generate circle
        x.append(x[0])
        y = [p[1] for p in line]
        y.append(y[0])
        if z is None:
            ax.plot(x, y, color=color, alpha=alpha, linewidth=linewidth, linestyle='-', label=label)
        else:
            ax.plot(x, y, [z] * (len(x)), color=color, alpha=alpha, linewidth=linewidth, linestyle='-', label=label)
        label = None
        G.remove_nodes_from(line)

    # draw existing edges
    exist_edges = list(G.edges())
    for edge in exist_edges:
        x = [edge[0][0], edge[1][0]]
        y = [edge[0][1], edge[1][1]]
        if z is None:
            ax.plot(x, y, color=color, alpha=alpha, linewidth=linewidth, linestyle='-', label=label)
        else:
            ax.plot(x, y, [z] * (len(x)), color=color, alpha=alpha, linewidth=linewidth, linestyle='-', label=label)
    # lines = []
    # for i in range(len(edf)):
    #     linex, liney = [], []
    #     stack = []
    #     if edfs.iloc[i, 'visited'] is False:
    #         edfs.iloc[i, 'visited'] = True
    #         st_vertex =  edfs.loc[i, ['x1', 'y1']]
    #         tar_vertex = edfs.loc[i, ['x2', 'y2']]
    #         linex.extend([st_vertex[0], tar_vertex[0]])
    #         liney.extend([st_vertex[1], tar_vertex[1]])
    #         selection1 = (edfs['x1'] == tar_vertex[0]) * (edfs['y1'] == tar_vertex[1])
    #         selection2 = (edfs['x2'] == tar_vertex[0]) * (edfs['y2'] == tar_vertex[1])
    #         next_edge1 = edfs[selection1]
    #         next_edge2 = edfs[selection2]
    #         next_edge1 = next_edge1[next_edge1['visited'] == False]
    #         next_edge2 = next_edge1[next_edge2['visited'] == False]
    #

