import networkx as nx
import math
try:
    import matplotlib.pyplot as plt
except:
    raise


def create2D_Grid_Graph(gridSize, show = False, showEdges = False, multipleRadiosPerNodeArg = False, numChannels = 1, linkCap = 1.0, commRangeArg = 1.0, infRange = 1.0):
    global pos, numRadiosPerNode
    G = nx.grid_2d_graph(gridSize, gridSize)
    #print "edges", G.edges()
    G = nx.convert.convert_to_directed(G)
    G = nx.convert_node_labels_to_integers(G, ordering="sorted", label_attribute='coor')
    if not showEdges:
        G.remove_edges_from(G.edges())
    posG = {}
    counter = 0
    for i in range(gridSize):
        for j in range(gridSize):
            posG[counter] = (i, j)
            counter += 1
    if(show):
        nx.draw(G, posG, node_size = 1000.0, font_size = 14, node_color = 'white')
        plt.savefig("/home/hmedal/Documents/2_msu/1_MSU_Projects/PAPER_JammingSpatialInterference/figures/3x3_noEdges.pdf") # save as jpg
        plt.show() # display
        
create2D_Grid_Graph(3, show = True)