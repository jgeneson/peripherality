import csv
import networkx as nx

def ng(G):
    """Computes n_G(u,v) for every ordered pair of vertices u,v of G.
    
    For every ordered pair of vertices u,v of graph G, computes n_G(u,v), the
    number of vertices of G that are closer to u than to v.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys (u,v) to n_G(u,v).
    """
    distance = dict(nx.all_pairs_shortest_path_length(G))
    ngraph = {}
    vertices = G.nodes
    for x in vertices:
        for y in vertices:
            ngraph[(x,y)] = 0
            for z in vertices:
                if distance[x][z] < distance[y][z]:
                    ngraph[(x,y)] += 1
    return ngraph
                        
                        
def peripherality(G):
    """Computes peri(u) for each vertex u of G.
    
    For every vertex u of graph G, computes peri(u), the number of vertices in
    G that has more vertices of G closer to itself than u.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys u to peri(u).
    """
    ngraph = ng(G)
    peri = {}
    vertices = G.nodes
    for x in vertices:
        peri[x] = 0
        for y in vertices:
            if ngraph[(y,x)] > ngraph[(x,y)]:
                peri[x] += 1
                    
    return peri

def spr(G):
    """Computes spr(u) for each vertex u of G.
    
    For every vertex u of graph G, computes spr(u), the sum of n_G(x,u) where
    x runs over all vertices of G.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys u to spr(u).
    """
    ngraph = ng(G)
    sprt = {}
    vertices = G.nodes
    for x in vertices:
        sprt[x] = 0
        for y in vertices:
            sprt[x] += ngraph[(y,x)]
                    
    return sprt

def eperi(G):
    """Computes eperi(e) for each edge e of G.
    
    For every edge e = (u,v) of graph G, computes eperi(e), the number of
    vertices x of G such that n_G(u,x) < n_G(x,u) and n_G(v,x) < n_G(x,v).

    Args:
        G: a graph.
    Returns:
        A dict mapping keys e to eperi(e).
    """
    ngraph = ng(G)
    eper = {}
    vertices = G.nodes
    edges = G.edges
    for e in edges:
        eper[e] = 0
        for x in vertices:
            if x != e[0] and x != e[1] and ngraph[(x,e[0])] > ngraph[(e[0],x)] and ngraph[(x,e[1])] > ngraph[(e[1],x)]:
                eper[e] += 1
            
    return eper

def eperi1(G):
    ngraph = ng(G)
    eper = {}
    vertices = G.nodes
    edges = G.edges
    for e in edges:
        eper[e] = 0
        for x in vertices:
            if x != e[0] and x != e[1] and (ngraph[(x,e[0])] > ngraph[(e[0],x)] or ngraph[(x,e[1])] > ngraph[(e[1],x)]):
                eper[e] += 1
            
    return eper

def espr(G):
    """Computes espr(e) for each edge e of G.
    
    For every edge e = (u,v) of graph G, computes espr(e), the sum of
    n_G(x,u) + n_G(x,v) with x running over all vertices of G except for u and
    v.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys e to espr(e).
    """
    ngraph = ng(G)
    esprt = {}
    vertices = G.nodes
    edges = G.edges
    for e in edges:
        esprt[e] = 0
        for x in vertices:
            if x != e[0] and x != e[1]:
                esprt[e] += (ngraph[(x,e[0])]+ngraph[(x,e[1])])
            
    return esprt

def irr(G):
    """Computes irr(e) for each edge e of G.
    
    For every edge e = (u,v) of graph G, computes irr(e), the difference of
    degrees of u and v.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys e to irr(e).
    """
    edges = G.edges
    degree = G.degree
    irre = {}
    for x in edges:
        irre[x] = abs(degree[x[0]]-degree[x[1]])
    return irre

def mostar(G):
    """Computes mostar(G) for graph G.
    
    For graph G, computes mostar(G), the sum of difference between n_G(u,v) and
    n_G(v,u) with (u,v) running over all edges of G.

    Args:
        G: a graph.
    Returns:
        mostar(G).
    """
    edges = G.edges
    ngraph = ng(G)
    most = {}
    for x in edges:
        most[x] = abs(ngraph[(x[0],x[1])]-ngraph[(x[1],x[0])])
    return most

def edeg(G):
    """Computes edeg(e) for every edge of graph G.
    
    For every edge e = (u,v) of graph G, computes edeg(G), the number of
    vertices of G not equal to u or v which are adjacent to u, v, or both.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys e to edge(e).
    """
    edges = G.edges
    vertices = G.nodes
    edegree = {}
    for x in edges:
        counter = 0
        for y in vertices:
            if y != x[0] and y != x[1] and ((x[0],y) in edges or (x[1],y) in edges):
                counter += 1
        edegree[x] = counter
    return edegree
            

def edge_ecc(G):
    """Computes edeg_ecc(e) for every edge of connected graph G.
    
    For every edge e = (u,v) of connected graph G, computes edeg_ecc(e), the
    longest distance from e to any vertex of G.

    Args:
        G: a graph.
    Returns:
        A dict mapping keys e to edeg_ecc(e).
    """
    edges = G.edges
    vertices = G.nodes
    distance = dict(nx.all_pairs_shortest_path_length(G))
    edge_ecc = {}
    for x in edges:
        maxdist = 0
        for y in vertices:
            maxdist = max(maxdist,min(distance[x[0]][y],distance[x[1]][y]))
        edge_ecc[x] = maxdist
    return edge_ecc
    

def listrank(index_list):
    """Computes the list rank of a dict.
    
    For every key x of a dict, computes the number of other keys with larger
    value.

    Args:
        index_list: a dict.
    Returns:
        A dict mapping a key x to the number of other keys with larger value in
	the input dict.
    """
    list_rank = {}
    for z in index_list:
        counter = 1
        for y in index_list:
            if index_list[z] < index_list[y]:
                counter += 1
        list_rank[z] = counter
    return list_rank


def revlistrank(index_list):
    """Computes the reverse list rank of a dict.
    
    For every key x of a dict, computes the number of other keys with smaller
    value.

    Args:
        index_list: a dict.
    Returns:
        A dict mapping a key x to the number of other keys with smaller value in
	the input dict.
    """
    list_rank = {}
    for z in index_list:
        counter = 1
        for y in index_list:
            if index_list[z] > index_list[y]:
                counter += 1
        list_rank[z] = counter
    return list_rank

reaction = {}
counter = 0

with open('downloads/super_fast_reactions.csv', mode = 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        counter += 1
        reactants = []
        products = []
        reactant_check = 0
        for j in range(len(row)):
            if row[j] != '.' and reactant_check == 0:
                reactants.append(row[j])
            elif row[j] == '.':
                reactant_check = 1
            elif row[j] != '.' and reactant_check == 1:
                products.append(row[j])
        reaction[counter] = [reactants, products]

reaction_graph_sf = nx.Graph()

for i in reaction:
    new_edge = []
    for j in reaction[i][0]:
        if j != 'M':
            reaction_graph_sf.add_node(j)
            new_edge.append(j)
    for s in range(len(new_edge)):
        for t in range(s+1,len(new_edge)):
            if new_edge[s] != new_edge[t]:
                reaction_graph_sf.add_edge(new_edge[s],new_edge[t])

reaction = {}
counter = 0

with open('downloads/mozart4.csv', mode = 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in reader:
        counter += 1
        if row[0] != row[2] and row[0] != 'M' and row[2] != 'M':
            reaction[counter] = [row[0],row[2]]
           
reaction_graph_m4 = nx.Graph()

for i in reaction:
    new_edge = []
    for j in reaction[i]:
        if j != 'M':
            reaction_graph_m4.add_node(j)
            new_edge.append(j)
    for s in range(len(new_edge)):
        for t in range(s+1,len(new_edge)):
            if new_edge[s] != new_edge[t]:
                reaction_graph_m4.add_edge(new_edge[s],new_edge[t])


# superfast calculations
def central_calc(G):
    peri_list = peripherality(G)

    spr_list = spr(G)

    dc = nx.degree_centrality(G)
    dc_rank = listrank(dc)

    cc = nx.closeness_centrality(G)
    cc_rank = listrank(cc)

    bc = nx.betweenness_centrality(G)
    bc_rank = listrank(bc)

    ec = nx.eigenvector_centrality(G)
    ec_rank = listrank(ec)

    ecc = nx.eccentricity(G)


    for z in sorted(peri_list):
        print z, '&', revlistrank(peri_list)[z], '&', revlistrank(spr_list)[z], '&', dc_rank[z], '&', cc_rank[z], '&', bc_rank[z], '&', ec_rank[z], '&', revlistrank(ecc)[z], ' \\\ '
        print '\hline'

    print(G.order())
    print(G.number_of_edges())
    print(nx.diameter(G))
    print(G.degree())

    eecc = edge_ecc(G)

    edegree = edeg(G)

    eperx = eperi(G)

    espr1 = espr(G)

    mst = mostar(G)



    erank = []
    for z in edegree:
        erankrow = []
        if z[0] < z[1]:
            erankrow.append(z[0])
            erankrow.append(z[1])
        else:
            erankrow.append(z[1])
            erankrow.append(z[0])
        erankrow.append(listrank(edegree)[z])
        erankrow.append(revlistrank(eecc)[z])

        erankrow.append(revlistrank(eperx)[z])
        erankrow.append(revlistrank(espr1)[z])
        erankrow.append(revlistrank(mst)[z])
        erank.append(erankrow)


    for z in sorted(erank):
        print z[0],',',z[1], '&', z[2],'&', z[3], '&', z[4], '&', z[5], '&', z[6], ' \\\ '
        print '\hline'
    print ''
    return ''

print(central_calc(reaction_graph_sf))
print(central_calc(reaction_graph_m4))
