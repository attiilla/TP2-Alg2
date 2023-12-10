import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import time
from datetime import datetime

filename = sys.argv[1][5:]
#print(sys.argv)

def distance(p0,p1):
    return np.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

def compute_value(g, h):
    v = 0
    for i in range(len(h)-1):
        v += g[h[i]][h[i+1]]['weight']
    v += g[h[0]][h[-1]]['weight']
    return v    

def estimate_state(g, prev_e=0, last_vertex = None, new_vertex = None):
    e = prev_e
    if last_vertex==None:
        for vid in g.nodes():
            sorted_edges = sorted(g.edges(vid, data=True), key=lambda x: x[2]['weight'])[:2]
            e += sorted_edges[0][2]["weight"]
            e += sorted_edges[1][2]["weight"]
        e = e/2
    else:
        #new
        new_edge = g.get_edge_data(new_vertex,last_vertex)
        sorted_edges = sorted(g.edges(new_vertex, data=True), key=lambda x: x[2]['weight'])[:2]
        if not (sorted_edges[0][2]["id"] == new_edge["id"]):
            e += (new_edge["weight"] - sorted_edges[1][2]["weight"])/2
        
        #last
        sorted_edges = sorted(g.edges(last_vertex, data=True), key=lambda x: x[2]['weight'])
        e += (g.get_edge_data(new_vertex,last_vertex)["weight"] - sorted_edges[0][2]["weight"])/2
    return e

def log_approximate(best,s,filename):
    log_filename = "./results/"+filename.split(".")[0]+"_log.txt"
    with open(log_filename, 'a') as log_file:
        log_file.write(f"{best[0]}\\{best[1]}\\{s}\n")

def log_tsp(best, max_expanded, filename, init_time):
    current_timestamp = time.time()
    current_datetime = datetime.fromtimestamp(current_timestamp)
    log_filename = "./results/"+filename.split(".")[0]+"_log.txt"
    with open(log_filename, 'r') as file:
        lines = file.readlines()
    if len(lines)==2:
        lines.append(f"{init_time.strftime('%Y-%m-%d %H:%M:%S.%f')}\\{current_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')}\\{best[0]}\\{best[1]}\\{max_expanded}\n")
    else:
        lines[2] = f"{init_time.strftime('%Y-%m-%d %H:%M:%S.%f')}\\{current_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')}\\{best[0]}\\{best[1]}\\{max_expanded}\n"

    with open(log_filename, 'w') as file:
        file.writelines(lines)

def build_graph_from_file(filename):
    location = "./Data/"+filename
    g = nx.Graph()
    file = open(location,'r')
    i = 0
    n = 0
    while True:
        line = file.readline()
        if i==3:
            n = int(line.replace(" ", "").split(":")[1])
            g.add_nodes_from(range(n))
        if i==5:
            break
        i+=1
    x = []
    y = []
    node_data = file.readline()
    while node_data!='EOF\n':
        arr = node_data[:-1].split()
        x.append(arr[1])
        y.append(arr[2])
        node_data = file.readline()
        eid = 0
    for i in range(n):
        for j in range(i+1,n):
            p0 = [float(x[i]), float(y[i])]
            p1 = [float(x[j]), float(y[j])]
            g.add_edge(i, j, weight=distance(p0, p1), id = eid)
            eid+=1
    return g

def TwiceAround(g, filename):
    T = nx.algorithms.tree.minimum_spanning_tree(g, algorithm='prim')
    h = list(nx.dfs_preorder_nodes(T))
    ans = (h, compute_value(g,h))
    log_approximate(ans,2*len(T.nodes())-1,filename)
    return ans

def Christofides(g, filename):
    T = nx.algorithms.tree.minimum_spanning_tree(g, algorithm='prim')
    TM = nx.MultiGraph(T)
    odds = []
    for node in T.nodes():
        if T.degree(node)%2==1:
            odds.append(node)
    H = g.subgraph(odds)
    matching = nx.algorithms.matching.min_weight_matching(H)
    for e in matching:
        TM.add_edge(e[0],e[1])
    h = list(nx.dfs_preorder_nodes(TM))
    ans = (h, compute_value(g,h))
    log_approximate(ans,4*len(T.nodes())-2+len(matching),filename)
    return ans

def recursive_branch_and_bound(g, vid, e, h, best, n, depth, filename, init_time, max_expanded = 0):
    depth+=1
    if depth>max_expanded:
        max_expanded = depth
    if len(h)==n:
        value = compute_value(g,h)
        if value<best[1]:
            best = (list(h),value)
        return best, max_expanded
    edges = sorted(g.edges(vid, data=True), key=lambda x: x[2]['weight'])
    for edg in edges:
        if not (edg[1] in h):
            e = estimate_state(g, prev_e = e, last_vertex = vid, new_vertex = edg[1])
            if e<best[1]:
                h.append(edg[1])
                b, max = recursive_branch_and_bound(g, edg[1], e, h, best, n, depth, filename, init_time, max_expanded=max_expanded)
                h.pop()
                if max>max_expanded:
                    max_expanded = max
                    log_tsp(best, max_expanded, filename, init_time)
                if b[1]<best[1]:
                    best = b
                    log_tsp(best, max_expanded, filename, init_time)
    return best, max_expanded

def branch_and_bound_tsp(g, filename, vid = 0):
    init_time = datetime.now()
    n = nx.number_of_nodes(g)
    cycle = list(range(n))
    best = (cycle, compute_value(g,cycle))
    h = [vid]
    e = estimate_state(g, prev_e = 0, last_vertex = None, new_vertex = None)
    depth = 1
    best, max_expanded = recursive_branch_and_bound(g, vid, e, h, best, n, depth, filename, init_time)
    log_tsp(best, max_expanded, filename, init_time)
    return best, max_expanded

G = build_graph_from_file(filename)
TwiceAround(G,filename)
Christofides(G,filename)
branch_and_bound_tsp(G,filename)