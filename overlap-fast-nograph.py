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

def get_merged_node_maximal_overlap(reads, step):                  
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



