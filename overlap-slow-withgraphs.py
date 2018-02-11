def overlap(a,b, min_length = 1):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if(start == -1):
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

from itertools import permutations
from graphviz import Digraph
import os

def get_merged_node_maximal_overlap(reads, step):
    dot = Digraph(format="pdf")
    
    graph = []
    for _ in reads:
        graph.append([])

    for i,a in enumerate(reads):
        for b in reads:
            graph[i].append(overlap(a,b))        
    
    for i,read in enumerate(reads):
        dot.node(read, read)
    for i,a in enumerate(reads):
        for j,b in enumerate(reads):
            if(i is not j and graph[i][j] is not 0):
                dot.edge(a,b, label=str(graph[i][j]))
                
    dot.render('overlap-output/step-' + str(step), cleanup=True)                  

    read_a = None; read_b = None;
    overlap_len = 0
    for a,b in permutations(reads, 2):
        o_len = overlap(a,b)
        if(o_len > overlap_len):
            overlap_len = o_len
            read_a, read_b = a,b
    return read_a, read_b, overlap_len

def overlap_reads(reads):
    reads = sorted(reads)
    step = 1
    if os.path.exists("overlap-output/"):
        os.system("rm -rf overlap-output/")
    reada, readb , overlap = get_merged_node_maximal_overlap(reads, step)
    while overlap > 0:
        reads.remove(reada)
        reads.remove(readb)        
        reads.append((reada + readb[overlap:]))
        step += 1
        reada, readb , overlap = get_merged_node_maximal_overlap(reads,step)
    return ''.join(reads)

a = int(input("Enter number of reads: "))
b = []
for x in range(a):
    b.append(input("Enter read" + str(x + 1) + " : "))

print("The obtained superstring is " + overlap_reads(b))
print("Overlap graphs generated stepwise are in overlap-output directory.")


