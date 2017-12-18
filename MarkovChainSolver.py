'''
@author J.D. Mazz
This program solves for the steady state of any arbitrary probability graph. 
This program can handle cyclic Digraphs with/without terminal nodes and 
determine the probability distributions of all nodes within the network.
Takes a matrix of counts, the desired level of convergence (necessary for cycles), 
and the max precision of the denominator of the fractional solution and returns
The counts indicate the ratios of the transtions and the column sums will become
the denominators of the transition matrix. Each row of the transition matrix
will sum to one. Once the steady state is achieved, the row vector solution 
will be normalized and indicates the relative frequency of each node.
Ex: [terminal_node_numerator1, ..., common_denominator]
'''
from fractions import Fraction

def MarkovChainSolver(m, epsilon, denominator_cap):
    T = make_trans_mat(m)
    v = T[0]
    oldv = [-1]*len(v)
    while not within(v, oldv, epsilon):
        oldv = v
        v = vec_x_mat(v, T)
    norm = sum(v)**-1
    v = [x*norm for x in v]
    soln = reformat(v, denominator_cap)
    terms = end_nodes(m)
    fsoln = []
    for i in terms:
        fsoln.append(soln[i])
    fsoln.append(soln[-1])
    return fsoln

def vec_x_mat(v, mat):
    newv = [0]*len(v)
    for r in range(len(v)):
        rsum = 0
        for c in range(len(v)):
            rsum += v[c]*mat[c][r]
        newv[r] = rsum
    return newv

def make_trans_mat(m):
    T = [[0]*len(m) for i in range(len(m))]
    for i in range(len(m)):
        rsum = sum(m[i])
        if rsum == 0:
            T[i] = [Fraction(0,1)]*len(m)
            T[i][i] = Fraction(1,1) # diag soln
            continue
        for j in range(len(m)):
            T[i][j] = Fraction(m[i][j],rsum)
    return T

def gcd(x, y):
    while y:      
        x, y = y, x % y
    return x

def lcm(x, y):
    return x * y // gcd(x, y)

def lcd(probs):
    if len(probs) == 0:
        return 0
    if len(probs) == 1:
        return probs[0].denominator
    d = lcm(probs[1].denominator, probs[0].denominator)
    for i in range(2, len(probs)):
        if i == 0:
            continue
        d = lcm(probs[i].denominator, d)
    return d

def reformat(p, denominator_cap):
    out = []
    p = [f.limit_denominator(denominator_cap) for f in p]
    d = lcd(p)
    for f in p:
        c = d // f.denominator
        out.append(c * f.numerator)
    out.append(d)
    return out

def end_nodes(m):
    terms = []
    for i in range(len(m)):
        if sum(m[i]) == 0:
            terms.append(i)
    return terms

def within(v, newv, e):
    diff = 0
    for i in range(len(v)):
        diff += abs(v[i] - newv[i])
    if diff > e:
        return False
    return True

def test():
    M = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    print(MarkovChainSolver(M, 0.001, 100))