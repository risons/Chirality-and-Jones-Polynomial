from sympy import symbols,simplify
from copy import deepcopy
A,t = symbols('A t',positive = True)   
def kauffman(pairs):
    crossings = [pair[0]for pair in pairs]
    print(pairs)
    if len(crossings) == 0 or (len(crossings) == 1 and len(crossings[0]) == 0):
        #trivial case
        return 1
    
    #working cross by cross
    cross = crossings[0]
    if sign(cross,crossings) != pairs[0][1]:
        print(pairs)
        print(crossings,pairs[0][1],sign(cross,crossings))
    #initiallising each element of cross
    a,b,c,d = cross
    connections = [0 if len(arcConnect(crossings,i)) ==1 else crossings.index(arcConnect(crossings,i)[0]) if arcConnect(crossings,i)[1] == cross else crossings.index(arcConnect(crossings,i)[1]) for i in cross]
    E,F,G,H = connections
    testA = aSmoothing(pairs,a,b,c,d,E,F,G,H)
    testB = bSmoothing(pairs,a,b,c,d,E,F,G,H)
    if E==0 and F == 0 and G == 0 and H== 0:
        # if both arcs are loops
        if len(crossings)==1:
            if pairs[0][1] == 1: # if the western and eastern arcs are loops        
                return -A **3 
            else: 
                return -A**-3    # if the northern and southern arcs are loops
        else:
            if pairs[0][1] == 1:
                return (A**5 + A)* kauffman(pairs[1:])
            else:
                return (A**-5 + A**-1)* kauffman(pairs[1:])  
    elif  (G==0 and F ==0) or (E == 0 and H == 0): # if the northern or southern arcs are loops           
        return simplify(kauffman(testB) * (-A**1 - A**-3) + A* kauffman(testA))
    elif (G==0 and H ==0) or (E == 0 and F == 0): # if the western or eastern arcs are loops 
        return simplify(kauffman(testA) * (-A**3 - A**-1) + (A**-1)* kauffman(testB))
    else:
        return simplify(A* kauffman(testA) + (A**-1)* kauffman(testB))
    
def aSmoothing(originalpairs,a,b,c,d,E,F,G,H):
    #connecting a (se) to b (ne) and d (sw) to c (nw)
    pairs = deepcopy(originalpairs)
    crossings = [pair[0] for pair in pairs]
    signs = [pair[1] for pair in pairs]
    if G == 0 and F == 0:
        h = crossings[H]
        h[h.index(d)] = a
    elif E == 0 and H == 0:
        f = crossings[F]
        f[f.index(b)] = c
    else:
        if E!= 0 or F != 0: 
            f = crossings[F]
            f[f.index(b)] = a
        if G != 0 or H !=0: 
            h = crossings[H]
            h[h.index(d)] = c
    del crossings[0]
    del signs[0]

    #relabeling
    arcs = sorted({x for crossing in crossings for x in crossing})
    mapping = {old:new+1 for new,old in enumerate(arcs)}
    crossings = [[mapping[x] for x in crossing] for crossing in crossings]
    return [(crossings[i],signs[i])for i in range(len(crossings))]

def bSmoothing(originalpairs,a,b,c,d,E,F,G,H):  
    pairs = deepcopy(originalpairs)
#connecting a (se) to d (sw) and b (ne) to c (nw)
    crossings = [pair[0] for pair in pairs]
    signs = [pair[1] for pair in pairs]
    if G == 0 and H == 0:
        f = crossings[F]
        f[f.index(b)] = a
    elif E == 0 and F == 0:
        g = crossings[G]
        g[g.index(c)] = d
    else:
        if G != 0 or F!= 0: 
            g = crossings[G]
            g[g.index(c)] = b        
        if H !=0 or E != 0: 
            h = crossings[H]

            h[h.index(d)] = a
 
    del crossings[0]
    del signs[0]
    #relabeling
    arcs = sorted({x for crossing in crossings for x in crossing})
    mapping = {old:new+1 for new,old in enumerate(arcs)}
    crossings = [[mapping[x] for x in crossing] for crossing in crossings]
    return [(crossings[i],signs[i])for i in range(len(crossings))]       
def arcConnect(crossings,arc):
    return [cross for cross in crossings if arc in cross]
def succ(arc,crossings):
    succs = {}
    for comp in get_components(crossings):
        for i in range(len(comp)):
            succs[comp[i-1]] = comp[i]
    return succs[arc]
def sign(crossing,crossings):
        if succ(crossing[3],crossings) == crossing[1]:
            if succ(crossing[1],crossings) == crossing[3]:
                if crossing[2] == crossing[3]:
                    return 1
                else:
                    return -1
            else:
                return 1
        
        else:
            return -1
def get_components(crossings):
        """Determines the components of a ),(
           Returns a list of lists of (adjacently ordered) arc labels """
    
        components = []
        # hash giving both neighbouring arc labels of an arc 
        nbrs = {}
        for c in crossings:
            nbrs.setdefault(c[0], []).append(c[2])
            nbrs.setdefault(c[2], []).append(c[0])
        for c in crossings:
            nbrs.setdefault(c[1], []).append(c[3])
            nbrs.setdefault(c[3], []).append(c[1])

        visited = set()
        # iterate over all arcs
        for arc in nbrs:
            if arc in visited:
                continue
            stack = [arc]
            visited.add(arc)
            knot = []
            
            # do depth-first-search on the arc to find its knot
            while stack:
                curr = stack.pop()
                knot.insert(0,curr)
                for nbr in nbrs[curr]:
                    if nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
            components.append(knot)

        return components
def Jones(pd):
    pairs = pd.pairs
    
    k = simplify((-A**-3)**(pd.writhe())* kauffman(pairs))
    if len(pd.get_components())==1:
        return simplify(k.subs(A**-4,t))
    else:
        return simplify(k.subs(A**-2,t))

class Diagram:
    def __init__(self, pd):
        # matrix is a 2d list of segment numbers
        self.pd = pd
        self.components = self.get_components()
        self.succs = {}
        self.arc_comps = {}
        for comp in self.components:
            for i in range(len(comp)):
                self.succs[comp[i-1]] = comp[i]
                self.arc_comps[comp[i]] = comp
        self.pairs = [(crossing,self.get_sign(crossing)) for crossing in self.pd]
        
    
    def get_sign(self,crossing):
        if self.get_arc_succ(crossing[3]) == crossing[1]:
            if self.get_arc_succ(crossing[1]) == crossing[3]:
                if crossing[2] == crossing[3]:
                    return 1
                else:
                    return -1
            else:
                return 1
        
        else:
            return -1

        
    def writhe(self):
        return sum([self.get_sign(c) for c in self.pd])

    def get_arc_comp(self,arc):
        return self.arc_comps[arc]
    
    def get_linking_number(self):
        acc = 0
        for crossing in self.pd:
            if self.get_arc_comp(crossing[0]) != self.get_arc_comp(crossing[1]):
                result += self.get_sign(crossing)
        return 1/2 * acc
                
                
    def get_arc_succ(self, arc):
        return self.succs[arc]
    
    
   
    
    
    def get_components(self):
        """Determines the components of a PD
           Returns a list of lists of (adjacently ordered) arc labels """
    
        components = []
        # hash giving both neighbouring arc labels of an arc 
        nbrs = {}
        for c in self.pd:
            nbrs.setdefault(c[0], []).append(c[2])
            nbrs.setdefault(c[2], []).append(c[0])
        for c in self.pd:
            nbrs.setdefault(c[1], []).append(c[3])
            nbrs.setdefault(c[3], []).append(c[1])

        visited = set()
        # iterate over all arcs
        for arc in nbrs:
            if arc in visited:
                continue
            stack = [arc]
            visited.add(arc)
            knot = []
            
            # do depth-first-search on the arc to find its knot
            while stack:
                curr = stack.pop()
                knot.insert(0,curr)
                for nbr in nbrs[curr]:
                    if nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
            components.append(tuple(knot))

        return components
print(kauffman([([1,5,2,4],1),([3,1,4,6],1),([5,3,6,2],1)]))
'''
links = [([ [4, 1, 3, 2],  [2, 3, 1, 4]],	-x **(-5)-x **(-1)
),([ [4, 2, 3, 1],  [2, 4, 1, 3]],	-x-x **5
),([ [6, 1, 7, 2],  [8, 3, 5, 4],  [2, 5, 3, 6],  [4, 7, 1, 8]],	-x **(-9)-x **(-5) + x **(-3)-x **(-1)
),([ [6, 2, 7, 1],  [8, 4, 5, 3],  [2, 8, 3, 7],  [4, 6, 1, 5]],	-x **3-x **7 + x **9-x **11
),([ [6, 1, 7, 2],  [10, 7, 5, 8],  [4, 5, 1, 6],  [2, 10, 3, 9],  [8, 4, 9, 3]],	x **(-7)-2/x **5 + x **(-3)-2/x + x-x **3
),([ [8, 2, 9, 1],  [10, 7, 5, 8],  [4, 10, 1, 9],  [2, 5, 3, 6],  [6, 3, 7, 4]],	x **(-7)-2/x **5 + x **(-3)-2/x + x-x **3
),([ [6, 1, 7, 2],  [10, 3, 11, 4],  [12, 8, 5, 7],  [8, 12, 9, 11],  [2, 5, 3, 6],  [4, 9, 1, 10]],	-x **(-9) + x **(-7)-3/x **5 + 2/x **3-2/x + 2*x-x **3
),([ [10, 2, 11, 1],  [6, 4, 7, 3],  [12, 10, 5, 9],  [8, 6, 9, 5],  [2, 12, 3, 11],  [4, 8, 1, 7]],	-x **3 + x **5-3*x **7 + 2*x **9-2*x **11 + 2*x **13-x **15
),([ [8, 1, 9, 2],  [12, 5, 7, 6],  [10, 3, 11, 4],  [4, 11, 5, 12],  [2, 7, 3, 8],  [6, 9, 1, 10]]	,-x **(-15) + x **(-13)-2/x **11 + 2/x **9-2/x **7 + x **(-5)-x **(-3)
),([ [10, 2, 11, 1],  [12, 6, 7, 5],  [8, 4, 9, 3],  [4, 8, 5, 7],  [2, 12, 3, 11],  [6, 10, 1, 9]],	-x **3 + x **5-2*x **7 + 2*x **9-2*x **11 + x **13-x **15
),([ [8, 1, 9, 2],  [2, 9, 3, 10],  [10, 3, 11, 4],  [12, 5, 7, 6],  [6, 7, 1, 8],  [4, 11, 5, 12]]	,-x **(-17) + x **(-15)-x **(-13) + x **(-11)-x **(-9)-x **(-5)
),([ [10, 2, 11, 1],  [2, 10, 3, 9],  [8, 4, 9, 3],  [12, 6, 7, 5],  [6, 12, 1, 11],  [4, 8, 5, 7]]	,-x + x **3-x **5 + x **7-x **9-x **13
),([ [6, 1, 7, 2],  [12, 8, 9, 7],  [4, 12, 1, 11],  [10, 5, 11, 6],  [8, 4, 5, 3],  [2, 9, 3, 10]]	,4-x **(-6) + 3/x **4-2/x **2-2*x **2 + 3*x **4-x **6
),([ [6, 2, 7, 1],  [12, 5, 9, 6],  [4, 12, 1, 11],  [10, 8, 11, 7],  [8, 3, 5, 4],  [2, 9, 3, 10]]	,4-x **(-6) + 3/x **4-2/x **2-2*x **2 + 3*x **4-x **6
),([ [6, 1, 7, 2],  [12, 7, 9, 8],  [4, 9, 1, 10],  [10, 6, 11, 5],  [8, 4, 5, 3],  [2, 12, 3, 11]]	,4-x **(-6) + 3/x **4-2/x **2-2*x **2 + 3*x **4-x **6
),([ [6, 2, 7, 1],  [12, 6, 9, 5],  [4, 9, 1, 10],  [10, 7, 11, 8],  [8, 3, 5, 4],  [2, 12, 3, 11]]	,4-x **(-6) + 3/x **4-2/x **2-2*x **2 + 3*x **4-x **6
),([ [6, 1, 7, 2],  [10, 3, 11, 4],  [12, 7, 9, 8],  [8, 11, 5, 12],  [2, 5, 3, 6],  [4, 9, 1, 10]]	,x **(-14)-x **(-12) + 3/x **10-x **(-8) + 3/x **6-2/x **4 + x **(-2)
),([ [6, 2, 7, 1],  [10, 3, 11, 4],  [12, 6, 9, 5],  [8, 12, 5, 11],  [2, 8, 3, 7],  [4, 9, 1, 10]]	,-1 + x **(-2) + 3*x **2-x **4 + 3*x **6-2*x **8 + x **10
),([ [6, 1, 7, 2],  [10, 4, 11, 3],  [12, 8, 9, 7],  [8, 10, 5, 9],  [2, 5, 3, 6],  [4, 12, 1, 11]]	,-1 + x **(-2) + 3*x **2-x **4 + 3*x **6-2*x **8 + x **10
),([ [6, 2, 7, 1],  [10, 4, 11, 3],  [12, 5, 9, 6],  [8, 9, 5, 10],  [2, 8, 3, 7],  [4, 12, 1, 11]]	,-1 + x **(-2) + 3*x **2-x **4 + 3*x **6-2*x **8 + x **10)
]
for (crossings,answer) in links:
    if Jones(Diagram(crossings))==simplify(answer):
        print(True)

#print(Jones([[2,10,3,9],[4,14,5,13],[6,12,7,11],[8,2,9,1],[10,8,11,7],[12,6,13,5],[14,4,1,3]]))
crossingss = [([[1,5,2,4],[3,1,4,6],[5,3,6,2]],	t+ t**3-t**4),([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]],	t**(-2)-t**(-1)+ 1-t+ t**2)
,([[2,8,3,7],[4,10,5,9],[6,2,7,1],[8,4,9,3],[10,6,1,5]],	t**2+ t**4-t**5+ t**6-t**7)
,([[1,5,2,4],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,7,10,6]],	t-t**2+ 2*t**3-t**4+ t**5-t**6)
,([[1,7,2,6],[3,10,4,11],[5,3,6,2],[7,1,8,12],[9,4,10,5],[11,9,12,8]],	t**(-2)-t**(-1)+ 2-2*t+ t**2-t**3+ t**4)
,([[1,8,2,9],[3,11,4,10],[5,1,6,12],[7,2,8,3],[9,7,10,6],[11,5,12,4]],	t**(-1)-1+ 2*t-2*t**2+ 2*t**3-2*t**4+ t**5)
,([[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]],	-t**(-3)+ 2*t**(-2)-2*t**(-1)+ 3-2*t+ 2*t**2-t**3)
,([[1,9,2,8],[3,11,4,10],[5,13,6,12],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]],	t**3+ t**5-t**6+ t**7-t**8+ t**9-t**10)
,([[2,10,3,9],[4,14,5,13],[6,12,7,11],[8,2,9,1],[10,8,11,7],[12,6,13,5],[14,4,1,3]],	t-t**2+ 2*t**3-2*t**4+ 2*t**5-t**6+ t**7-t**8)
,([[1,9,2,8],[3,11,4,10],[5,1,6,14],[7,13,8,12],[9,3,10,2],[11,5,12,4],[13,7,14,6]],	t**2-t**3+ 2*t**4-2*t**5+ 3*t**6-2*t**7+ t**8-t**9)
,([[2,10,3,9],[4,12,5,11],[6,14,7,13],[8,4,9,3],[10,2,11,1],[12,8,13,7],[14,6,1,5]],	t-2*t**2+ 3*t**3-2*t**4+ 3*t**5-2*t**6+ t**7-t**8)
,([[2,10,3,9],[4,2,5,1],[6,14,7,13],[8,12,9,11],[10,4,11,3],[12,6,13,5],[14,8,1,7]],	t**2-t**3+ 3*t**4-3*t**5+ 3*t**6-3*t**7+ 2*t**8-t**9)
,([[1,13,2,12],[3,9,4,8],[5,1,6,14],[7,10,8,11],[9,3,10,2],[11,6,12,7],[13,5,14,4]],	t**(-1)-2+ 3*t-3*t**2+ 4*t**3-3*t**4+ 2*t**5-t**6)
,([[1,10,2,11],[3,13,4,12],[5,14,6,1],[7,5,8,4],[9,2,10,3],[11,9,12,8],[13,6,14,7]],	t**(-4)-2*t**(-3)+ 3*t**(-2)-4*t**(-1)+ 4-3*t+ 3*t**2-t**3)
,([[1,9,2,8],[3,7,4,6],[5,12,6,13],[7,3,8,2],[9,1,10,16],[11,15,12,14],[13,4,14,5],[15,11,16,10]],	t**(-2)-t**(-1)+ 2-2*t+ 2*t**2-2*t**3+ t**4-t**5+ t**6)
,([[1,10,2,11],[3,13,4,12],[5,15,6,14],[7,1,8,16],[9,2,10,3],[11,9,12,8],[13,5,14,4],[15,7,16,6]],	1-t+ 2*t**2-2*t**3+ 3*t**4-3*t**5+ 2*t**6-2*t**7+ t**8)
,([[6,2,7,1],[14,10,15,9],[10,5,11,6],[12,3,13,4],[4,11,5,12],[2,13,3,14],[16,8,1,7],[8,16,9,15]],	t**(-4)-t**(-3)+ 2*t**(-2)-3*t**(-1)+ 3-3*t+ 2*t**2-t**3+ t**4)
,([[2,11,3,12],[4,9,5,10],[6,16,7,15],[8,14,9,13],[10,1,11,2],[12,3,13,4],[14,8,15,7],[16,6,1,5]],	t**(-5)-2*t**(-4)+ 3*t**(-3)-3*t**(-2)+ 3*t**(-1)-3+ 2*t**1-t**2+ t**3)
,([[1,7,2,6],[3,9,4,8],[5,12,6,13],[7,3,8,2],[9,15,10,14],[11,1,12,16],[13,4,14,5],[15,11,16,10]],	1-t+ 3*t**2-3*t**3+ 3*t**4-4*t**5+ 3*t**6-2*t**7+ t**8)
,([[2,9,3,10],[4,14,5,13],[6,16,7,15],[8,12,9,11],[10,1,11,2],[12,8,13,7],[14,6,15,5],[16,4,1,3]],	t**(-1)-1+ 3*t-4*t**2+ 4*t**3-4*t**4+ 3*t**5-2*t**6+ t**7)
,([[2,9,3,10],[4,14,5,13],[6,15,7,16],[8,1,9,2],[10,5,11,6],[12,4,13,3],[14,12,15,11],[16,7,1,8]],	-t**(-6)+ 2*t**(-5)-3*t**(-4)+ 4*t**(-3)-4*t**(-2)+ 4*t**(-1)-2+ 2*t-t**2)
,([[1,7,2,6],[3,12,4,13],[5,9,6,8],[7,3,8,2],[9,16,10,1],[11,14,12,15],[13,4,14,5],[15,10,16,11]],	-t**(-5)+ 2*t**(-4)-3*t**(-3)+ 4*t**(-2)-4*t**(-1)+ 5-3*t+ 2*t**2-t**3)
,([[6,2,7,1],[14,8,15,7],[10,3,11,4],[2,13,3,14],[12,5,13,6],[4,11,5,12],[16,10,1,9],[8,16,9,15]],	t**(-4)-2*t**(-3)+ 3*t**(-2)-4*t**(-1)+ 5-4*t+ 3*t**2-2*t**3+ t**4)
,([[2,14,3,13],[4,9,5,10],[6,11,7,12],[8,15,9,16],[10,5,11,6],[12,2,13,1],[14,7,15,8],[16,4,1,3]],	-t**(-6)+ 2*t**(-5)-4*t**(-4)+ 5*t**(-3)-4*t**(-2)+ 5*t**(-1)-3+ 2*t-t**2)
,([[1,10,2,11],[3,13,4,12],[5,15,6,14],[7,1,8,16],[9,2,10,3],[11,9,12,8],[13,7,14,6],[15,5,16,4]],	t**(-1)-2+ 4*t-4*t**2+ 5*t**3-5*t**4+ 3*t**5-2*t**6+ t**7)
,([[4,2,5,1],[10,8,11,7],[8,3,9,4],[2,9,3,10],[14,6,15,5],[16,11,1,12],[12,15,13,16],[6,14,7,13]],	t**(-4)-2*t**(-3)+ 4*t**(-2)-5*t**(-1)+ 5-5*t+ 4*t**2-2*t**3+ t**4)
,(	[[1,9,2,8],[3,14,4,15],[5,12,6,13],[7,11,8,10],[9,3,10,2],[11,16,12,1],[13,4,14,5],[15,6,16,7]],	-t**(-5)+ 2*t**(-4)-3*t**(-3)+ 5*t**(-2)-5*t**(-1)+ 5-4*t+ 3*t**2-t**3)
,	([[2,12,3,11],[4,8,5,7],[6,15,7,16],[8,14,9,13],[10,2,11,1],[12,10,13,9],[14,4,15,3],[16,5,1,6]],	t**(-1)-2+ 4*t-5*t**2+ 6*t**3-5*t**4+ 4*t**5-3*t**6+ t**7)
,	([[1,7,2,6],[3,15,4,14],[5,9,6,8],[7,3,8,2],[9,13,10,12],[11,1,12,16],[13,5,14,4],[15,11,16,10]],	t**2-2*t**3+ 5*t**4-5*t**5+ 6*t**6-6*t**7+ 4*t**8-3*t**9+ t**10)
,	([[2,7,3,8],[4,10,5,9],[6,1,7,2],[8,14,9,13],[10,15,11,16],[12,6,13,5],[14,3,15,4],[16,11,1,12]],	-t**(-6)+ 3*t**(-5)-5*t**(-4)+ 6*t**(-3)-6*t**(-2)+ 6*t**(-1)-4+ 3*t-t**2)
,	([[6,2,7,1],[14,8,15,7],[8,3,9,4],[2,13,3,14],[12,5,13,6],[4,9,5,10],[16,12,1,11],[10,16,11,15]],	t**(-4)-3*t**(-3)+ 5*t**(-2)-6*t**(-1)+ 7-6*t+ 5*t**2-3*t**3+ t**4)
,	([[6,2,7,1],[8,3,9,4],[16,11,1,12],[2,14,3,13],[4,15,5,16],[10,6,11,5],[12,7,13,8],[14,10,15,9]],	t**(-4)-4*t**(-3)+ 6*t**(-2)-7*t**(-1)+ 9-7*t+ 6*t**2-4*t**3+ t**4)
,	([[2,14,3,13],[5,11,6,10],[7,15,8,14],[9,5,10,4],[11,7,12,6],[12,2,13,1],[15,9,16,8],[16,4,1,3]],	t**3+ t**5-t**8)
,	([[1,7,2,6],[4,13,5,14],[5,9,6,8],[7,3,8,2],[10,15,11,16],[12,9,13,10],[14,3,15,4],[16,11,1,12]],	-t**(-5)+ t**(-4)-t**(-3)+ 2*t**(-2)-t**(-1)+ 2-t)
,	([[1,7,2,6],[4,13,5,14],[5,9,6,8],[7,3,8,2],[9,13,10,12],[11,1,12,16],[14,3,15,4],[15,11,16,10]],	2*t-2*t**2+ 3*t**3-3*t**4+ 2*t**5-2*t**6+ t**7)]
for (crossings,answer) in crossingss:
    if Jones(Diagram(crossings))==simplify(answer):
        print(True)
        
    else:
        print(False)
        print(Jones(Diagram(crossings)),"hi",simplify(answer))
        print(kauffman(crossings,0,Diagram(crossings).signs))
def links(crossings):
    linkList = []
    arcs = sorted({x for crossing in crossings for x in crossing})
    visited = set()
    for arc in arcs:
        if arc in visited:
            continue
        i = arc
        link = []
        while i not in link: 
            visited.add(i)
            link.append(i)
            for cross in crossings:
                if i in cross:
                    a,b,c,d = cross
                    if i ==a:
                        iTemp = c
                    elif i ==b:
                        iTemp = d
                    elif i ==c:
                        iTemp = a
                    else:
                        iTemp = b
                    if iTemp not in link:
                        i = iTemp
                        break
        linkList.append(link)
    return(linkList)
def crossingSign(cross,crossings):
    #assuming a is se, b is ne, c is nw and d is sw.
    a,b,c,d = cross
    linkList = links(crossings)
    bLink = [link for link in linkList if b in link][0]
    
    iLink = {i:[link for link in linkList if i in link][0] for i in cross}
    if any([iLink[i].index(i) == len(iLink[i])-1 and iLink[cross[cross.index(i)-2]].index(cross[cross.index(i)-2]) == 0 for i in cross]):
    #if bLink.index(b) == 0:
        if (b-d)*(c-a)>0:
            return -1
        else:
            return 1
    else:
        if (b-d)*(c-a)>0:
            return 1
        else:
            return -1

for (crossings,answer) in crossingss:
    arcs = sorted({x for crossing in crossings for x in crossing})
    print(all([Jones([[{old:(old+1)% len(arcs) if (old+1)% len(arcs) != 0 else len(arcs) for old in arcs}[x] for x in crossing] for crossing in crossings]) ==simplify(answer) for i in range(len(arcs))]))


crossings = [[3, 9, 4, 8], [5, 1, 6, 10], [7, 2, 1, 3], [2, 7, 8, 6], [9, 5, 10, 4]] #[[1,8,2,9],[3,11,4,10],[5,1,6,12],[7,2,8,3],[9,7,10,6],[11,5,12,4]]	
print(Jones(crossings))
print(writhe(crossings))
print([crossingSign(cross,crossings) for cross in crossings])

crossings = [[4, 1, 5, 8], [6, 2, 1, 3], [2, 6, 3, 5], [7, 4, 8, 7]] #[[2, 9, 3, 8], [4, 1, 5, 10], [6, 7, 7, 2], [1, 6, 8, 5], [9, 4, 10, 3]] #[[1,8,2,9],[3,11,4,10],[5,1,6,12],[7,2,8,3],[9,7,10,6],[11,5,12,4]]	
print(Jones(crossings))
print(writhe(crossings))
print([crossingSign(cross,crossings) for cross in crossings])

print([crossingSign(cross,crossings) for cross in crossings])
print(Jones(crossings))'''