from sage.all import *
from ore_algebra import OreAlgebra
from ore_algebra import nullspace
from ore_algebra.ore_algebra import OreAlgebra_generic

original_kronecker = nullspace.kronecker

# ore_algebra 0.5.0 has a Sage 10.9 compatibility typo: it constructs the
# associated multivariate commutative ring with one generator regardless of
# the number of Ore generators.  Keep the workaround local to this script.
def associated_commutative_algebra_fixed(self):
    try:
        return self._commutative_ring
    except AttributeError:
        self._commutative_ring = PolynomialRing(
            self.base_ring(), self.variable_names())
        return self._commutative_ring

OreAlgebra_generic.associated_commutative_algebra = \
    associated_commutative_algebra_fixed

def fraction_field_solver(mat, degrees=[], infolevel=0):
    """Kernel over Frac(base_ring), cleared back to the polynomial ring."""
    R0 = mat.base_ring()
    K0 = R0.fraction_field()
    ker = matrix(K0, mat).right_kernel_matrix()
    out = []
    for row in ker.rows():
        # The product is deliberately used instead of a multivariate lcm:
        # it is slower-growing only by content, and coercion is reliable.
        den = prod(a.denominator() for a in row)
        out.append(vector(R0, [R0(a * den) for a in row]))
    return nullspace._normalize(out)

def native_solver_fixed(self, R=None):
    return fraction_field_solver

OreAlgebra_generic._solver = native_solver_fixed

def kronecker_native_fixed(*args, **kwargs):
    # ore_algebra's fast Kronecker solver accepts ZZ[x,...] but rejects
    # QQ[x,...].  Clear the (constant) rational denominators first.
    fast = original_kronecker(nullspace.gauss())
    def solve(mat, degrees=[], infolevel=0):
        R0 = mat.base_ring()
        if R0.base_ring() is not QQ:
            return fast(mat, degrees=degrees, infolevel=infolevel)
        den = lcm([c.denominator() for a in mat.list()
                   for c in a.coefficients()] + [1])
        Rz = PolynomialRing(ZZ, R0.variable_names())
        mz = matrix(Rz, mat.nrows(), mat.ncols(),
                    [Rz(den*a) for a in mat.list()])
        ker = fast(mz, degrees=degrees, infolevel=infolevel)
        return [vector(R0, [R0(a) for a in row]) for row in ker]
    return solve

def quick_check_fixed(subsolver, *args, **kwargs):
    return subsolver

nullspace.kronecker = kronecker_native_fixed
nullspace.quick_check = quick_check_fixed

R = PolynomialRing(QQ, names=("n", "p", "q", "v"))
n,p,q,v = R.gens()
K = R.fraction_field()
A = OreAlgebra(K, "Sn", "Dp", "Dq", "Dv")
Sn,Dp,Dq,Dv = A.gens()

D = p*q*(1+v^2)+2*v
S = p^2*q^2*(1-p^2)*(1-q^2)*v^2

logp = (2*n+6)/p - 2*n*p/(1-p^2) - (2*n+4)*D.derivative(p)/D
logq = (2*n+5)/q - 2*(n+1)*q/(1-q^2) - (2*n+4)*D.derivative(q)/D
logv = (2*n+3)/v - (2*n+4)*D.derivative(v)/D

ideal = A.ideal([Sn-S/D^2, Dp-logp, Dq-logq, Dv-logv])
print("dimension",ideal.vector_space_dimension(),flush=True)

tv = ideal.ct(Dv, certificates=False, infolevel=1)
print("after v",len(tv),flush=True)
for op in tv: print(op,flush=True)
Av = tv[0].parent()
iv = Av.ideal(tv)
gv = dict(zip(Av.variable_names(), Av.gens()))

tq = iv.ct(gv["Dp"], certificates=False, infolevel=1)
print("after p",len(tq),flush=True)
for op in tq: print(op,flush=True)
Aq = tq[0].parent()
iq = Aq.ideal(tq)
gq = dict(zip(Aq.variable_names(), Aq.gens()))

tp = iq.ct(gq["Dq"], certificates=False, infolevel=2,
           early_termination=True)
print("after q",len(tp),flush=True)
for op in tp: print(op,flush=True)
