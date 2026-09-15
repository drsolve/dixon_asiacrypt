"""Exact checks for the worked examples. Python standard library only."""
from itertools import permutations, product
from fractions import Fraction

N = 5  # x, y, u, v, z
class P(dict):
    def __add__(a,b):
        b=poly(b); out=P(a)
        for m,c in b.items(): out[m]=out.get(m,0)+c
        return P({m:c for m,c in out.items() if c})
    __radd__=__add__
    def __neg__(a): return P({m:-c for m,c in a.items()})
    def __sub__(a,b): return a+-poly(b)
    def __rsub__(a,b): return poly(b)+-a
    def __mul__(a,b):
        out=P()
        for m,c in a.items():
            for n,d in poly(b).items():
                k=tuple(i+j for i,j in zip(m,n));out[k]=out.get(k,0)+c*d
        return P({m:c for m,c in out.items() if c})
    __rmul__=__mul__
    def __pow__(a,n):
        out=poly(1)
        for _ in range(n):out=out*a
        return out
    def at(a,values):return sum(c*prod(values[i]**m[i] for i in range(N)) for m,c in a.items())
def prod(xs):
    r=1
    for x in xs:r*=x
    return r
def poly(x):return x if isinstance(x,P) else P({(0,)*N:x}) if x else P()
def var(i):return P({tuple(int(j==i) for j in range(N)):1})
def det(A):
    n=len(A);r=poly(0)
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        r+=sign*prod(A[i][p[i]] for i in range(n))
    return r
x,y,u,v,z=map(var,range(N))
f=[x+y-z,x-y-1,x*x+y*y-5]
C=[[x+y-z,x-y-1,x*x+y*y-5],[u+y-z,u-y-1,u*u+y*y-5],[u+v-z,u-v-1,u*u+v*v-5]]
Chat=[f,[-1,-1,-(x+u)],[-1,1,-(y+v)]]
Delta=10-(z+1)*x+(1-z)*y+(2*x-z-1)*u+(2*y-z+1)*v
M=[[10,-z-1,1-z],[-z-1,2,0],[1-z,0,2]]
assert det(C)==(x-u)*(y-v)*Delta
assert det(Chat)==Delta
assert sum([1,x,y][i]*M[i][j]*[1,u,v][j] for i in range(3) for j in range(3))==Delta
assert det(M)==-4*(z*z-9)
for sol in [(2,1,3),(-1,-2,-3)]:
    vals=[sol[0],sol[1],0,0,sol[2]]
    assert all(g.at(vals)==0 for g in f)
    assert all(sum([1,sol[0],sol[1]][i]*poly(M[i][j]).at(vals) for i in range(3))==0 for j in range(3))
solutions=[(a,b,c) for a,b,c in product(range(17),repeat=3) if (a+b-c)%17==0 and (a-b-1)%17==0 and (a*a+b*b-5)%17==0]
assert solutions==[(2,1,3),(16,15,14)]
# One-variable Dixon counterexample: f1=t*x-1, f2=t*x-2; use z as t.
C2=[[z*x-1,z*x-2],[z*u-1,z*u-2]]
assert det(C2)==(x-u)*(-z)
print('PASS: cancellation identity, divided differences, coefficient matrix, determinant, kernel vectors, rational roots, all F_17 solutions, and false-converse example.')
