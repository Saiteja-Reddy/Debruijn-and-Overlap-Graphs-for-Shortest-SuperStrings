from graphviz import Digraph

def overlap(a,b, min_length = 1):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if(start == -1):
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def get_kmers(reads,k):
    kmers = set()
    edges = []
    for read in reads:
        for i in range(len(read) - k + 1):
            edges.append((read[i:i+k-1],read[i+1:i+k]))
            kmers.add(read[i:i+k-1])
            kmers.add(read[i+1:i+k])
    return kmers, edges

def get_debruijn_graph(reads,k):
    kmers, edges = get_kmers(reads,k)
    dot = Digraph(format="pdf")
    l_kmers = list(kmers)

    for i,read in enumerate(l_kmers):
        dot.node(read, read)
 
    for (a,b) in edges:
    	dot.edge(a,b)
    dot.render('debruijn-output/debruijn' , cleanup=True)
    return l_kmers, edges

def run_path(reads, k):
	kmers, edges =  get_debruijn_graph(["a_long_long_long", "ng_long_l", "g_long_time"],8)

	graphf = []
	graphb = []
	for _ in kmers:
	    graphf.append([])
	    graphb.append([])
	    
	for i,a in enumerate(kmers):
	    for b in kmers:
	        graphf[i].append(0)    
	        graphb[i].append(0)            
	        
	for (a,b) in edges:
	    a_i = kmers.index(a)
	    b_i = kmers.index(b)
	    graphf[a_i][b_i] += 1
	    graphb[b_i][a_i] += 1
	    
	starti = -1    
	for i in range(len(kmers)):
	    flag = 0
	    for j in range(len(kmers)):
	        if(graphb[i][j] > 0):
	            flag = 1
	    if(flag == 0):
	        starti = i
	        break
	        
	print("Starting walk from node: " + kmers[starti])


	joinmers = []
	flag = 1
	now = starti
	while flag:
	    flag=0
	    joinmers.append(kmers[now])
	    for j in range(len(kmers)):
	        if(graphf[now][j] > 0):
	            graphf[now][j] -= 1
	            flag = 1
	            now = j
	            break

	out = []
	out.append(joinmers[0])
	for i in range(1,len(joinmers)):
	    olen = overlap(joinmers[i-1],joinmers[i])
	    out.append(joinmers[i][olen:])

	print("The obtained superstring is " + ''.join(out))
	print("The Debruijn graph generated is in debruijn-output directory.")


a = int(input("Enter number of reads: "))
b = []
for x in range(a):
    b.append(input("Enter read" + str(x + 1) + " : "))
k = int(input("Enter k-mer size: "))

run_path(b,k)
