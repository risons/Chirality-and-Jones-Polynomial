from sympy import symbols,simplify
import timeit
A,t,x = symbols('A t x',positive = True)
def kauffman(crossings):
    if crossings is None or (len(crossings) == 1 and len(crossings[0]) == 0):
        #trivial case
        return 1
    
    #working cross by cross
    cross = crossings[0]
    #initiallising each element of cross
    a,b,c,d = cross
    #E, F, G, H store the indices of the crossings that arcs a, b, c, d connect to, respectively.
    #If an arc forms a loop (connects back to the same crossing), the index is 0.
    connections = [0 if len(arcConnect(crossings,i)) ==1 else crossings.index(arcConnect(crossings,i)[0]) if arcConnect(crossings,i)[1] == cross else crossings.index(arcConnect(crossings,i)[1]) for i in cross]
    E,F,G,H = connections 

    smoothedA = aSmoothing(crossings,a,b,c,d,E,F,G,H)
    smoothedB = bSmoothing(crossings,a,b,c,d,E,F,G,H)
    if E==0 and F == 0 and G == 0 and H== 0:
        # if both arcs are loops
        if len(crossings)==1: #if the only component of the diagram is this
            if cross[0] == cross[1]: # if the western and eastern arcs are loops        
                return -A **3 
            else: 
                return -A**-3    # if the northern and southern arcs are loops
        else: 
            if cross[0] == cross[1]:
                return (A**5 + A)* kauffman(crossings[1:])
            else:
                return (A**-5 + A**-1)* kauffman(crossings[1:])  
    elif  (G==0 and F ==0) or (E == 0 and H == 0): # if the northern or southern arcs are loops and multiplies that term by the necessary (-A^2 -A^-2) factor          
        return (kauffman(smoothedB) * (-A**1 - A**-3) + A* kauffman(smoothedA) )
    elif (G==0 and H ==0) or (E == 0 and F == 0): # if the western or eastern arcs are loops and multiplies that term by the necessary (-A^2 -A^-2) factor     
        return (kauffman(smoothedA) * (-A**3 - A**-1) + (A**-1)* kauffman(smoothedB))
    else:
        return (A* kauffman(smoothedA) + (A**-1)* kauffman(smoothedB))
    
def aSmoothing(originalcrossings,a,b,c,d,E,F,G,H):
    #Returns a list of crossings where the crossing [a,b,c,d] is smoothed such that connecting a (se) to b (ne) and d (sw) to c (nw)

    crossings = [crossing[:] for crossing in originalcrossings] #deep copy of originalcrossings

    if G == 0 and F == 0: #Checks if there a north loop and smooths a (se) goes to d (sw)
        h = crossings[H]
        h[h.index(d)] = a
    elif E == 0 and H == 0: #Checks if there isnt a south loop and smooths b (ne) to c (nw)
        f = crossings[F]
        f[f.index(b)] = c
    else:
        if E!= 0 or F != 0: #Checks if there is a east loop and smooths a (se) goes to b (ne)
            f = crossings[F]
            f[f.index(b)] = a
        if G != 0 or H !=0:  #Checks if there is a west loop and smooths d (sw) goes to c (nw) 
            h = crossings[H]
            h[h.index(d)] = c
    del crossings[0]
    '''
    #relabeling
    arcs = sorted({arc for crossing in crossings for arc in crossing})
    mapping = {old:new+1 for new,old in enumerate(arcs)}
    crossings = [[mapping[arc] for arc in crossing] for crossing in crossings]'''
    return crossings

def bSmoothing(originalcrossings,a,b,c,d,E,F,G,H):  
#Returns a list of crossings where the crossing [a,b,c,d] is smoothed such that a (se) goes to d (sw) and b (ne) to c (nw)

    crossings = [crossing[:] for crossing in originalcrossings] #deep copy of originalcrossings

    if G == 0 and H == 0: #Checks if there is a west loop and smooths a (se) goes to b (ne)
        f = crossings[F]
        f[f.index(b)] = a
    elif E == 0 and F == 0: #Checks if there is a east loop and smooths d (sw) goes to c (nw)
        g = crossings[G]
        g[g.index(c)] = d
    else:
        if G != 0 or F!= 0: #Checks if there isnt a north loop and smooths b (ne) to c (nw)
            g = crossings[G]
            g[g.index(c)] = b        
        if H !=0 or E != 0: #Checks if there isnt a south loop and smooths a (se) goes to d (sw)
            h = crossings[H]
            h[h.index(d)] = a
 
    del crossings[0]
    '''
    #relabeling
    arcs = sorted({arc for crossing in crossings for arc in crossing})
    mapping = {old:new+1 for new,old in enumerate(arcs)}
    crossings = [[mapping[arc] for arc in crossing] for crossing in crossings]'''
    return crossings      
def arcConnect(crossings,arc):
    # returns which crossing(s) a specific arc appears in 
    return [cross for cross in crossings if arc in cross]
def jones(pd):
    #Determines Jones polynomial from a planar diagram. Returns in terms of t.
    crossings = pd.pd
    k = ((-A**-3)**(pd.writhe())* kauffman(crossings))
    if len(pd.components)==1: # gives t = A^-4 for knots
        return simplify(k.subs(A**-4,t))
    else:
        return simplify(k.subs(A**-2,x)) # gives t = A^-2 for links
class PlanarDiagram:
    def __init__(self, pd):
        """pd is a list of crossings 
        A crossing is a list of 4 arc labels ordered anti-clockwise
        assume orientation of 1st -> 3rd arc of crossing being incoming underpass"""
        self.pd = pd
        self.components = self.find_components()
        
        self._succs = {} # hash giving next arc in an oriented strand
        self._arc_comps = {} # hash giving the component an arc belongs to
        
        for comp in self.components:
            for i in range(len(comp)):
                self._succs[comp[i-1]] = comp[i]
                self._arc_comps[comp[i]] = comp
        
    
    def get_sign(self,crossing):
        # finds the crossing sign based on succession of arcs in overpass
        if self.get_arc_succ(crossing[3]) == crossing[1]:
            if self.get_arc_succ(crossing[1]) == crossing[3]:
                # special 1 component knot case
                if crossing[2] == crossing[3]:
                    return 1
                else:
                    return -1
            else:
                return 1
        
        else:
            return -1

        
    def writhe(self):
        # sums signs of all crossings
        return sum([self.get_sign(c) for c in self.pd])


    def get_linking_number(self):
        # sums signs of all crossings that are made up of two components
        acc = 0
        for crossing in self.pd:
            if self.get_arc_comp(crossing[0]) != self.get_arc_comp(crossing[1]):
                result += self.get_arc_sign(crossing)
        return 1/2 * acc
                
                
    def get_arc_succ(self, arc):
        return self._succs[arc]
    
    
    def get_arc_comp(self,arc):
        return self._arc_comps[arc]

    
    def find_components(self):
        """Determines the components of a PD
           gives a list of tuples of (consecutively ordered) arc labels """
    
        components = []
        # hash giving consecutive arc labels of an arc in a strand in both directions
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
'''
for crossings, answer in crossingss:
    timer = timeit.Timer(
        lambda: jones(PlanarDiagram(crossings))
    )
    times = timer.repeat(repeat=3, number=1)
    print(times)
'''
for (crossings,answer) in crossingss:
    if simplify(answer.subs(t**-1,t)) == answer:
        print(True)
    else:
        print(False)
