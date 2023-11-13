from kummer_isogenies.gluing_isogeny_kummer import *
from kummer_isogenies.formulae_kummer import *
from kummer_isogenies.doubling_kummer import *
from sage.all import EllipticCurve



def compute_delta(type1_invariants,kernel2):
    #we assume that kernel = <((x-alpha),0), x*(x-beta, 0)>,
    #and then the third factor is (x-1/alpha)*(x-gamma)
    [A,B,C,E] = type1_invariants
    alpha_inv = kernel2[0][1]/kernel2[0][2]
    beta = kernel2[1][1]/kernel2[1][0]
    gamma = B - beta
    M = Matrix([kernel2[0][:3], kernel2[1][:3], [1,gamma+alpha_inv, gamma*alpha_inv]])
    return M.det()

def splitting_kummer_isogeny(ker4, phi):
    """
    input is a 4-torsion kernel lying above the next 2-isogeny, 
    and data of the previous (2,2)-isogeny
    We assume that the next (2,2)-kernel is of the form <((x-alpha),0), x*(x-beta, 0)>
    """
    invariants = phi[0]
    type1_invariants = invariants[-1]
    [Tp1,Tp2] = ker4
    new_Tp1 = double(Tp1,1,phi)
    new_Tp2 = double(Tp2,1,phi)
    ker2 = [new_Tp1,new_Tp2]
    [new_type1_invariants, special_trafo] = find_special_trafo(type1_invariants, ker2, Tp1)
    trafo_12 = compute_trafo_12(new_type1_invariants) #this is for transformation from type 1 to 2
    [A,B,C,E] = new_type1_invariants
    assert C == 1, "the isogeny is not split"

    #splitting maps as in the paper (Prop. 4.2), we ignore twists since we cannot recover the y-coordinate anyways
    #W-points of the elliptic curves are alpha1,alpha2,1 
    # and beta1,beta2,1 respectively
    AA = A*A
    BB = B*B
    denom = 1/((AA - 4)*(BB - 4))
    alpha1 = (AA + 4*A + 4)*(BB - 4)*denom
    alpha2 = (BB + 4*B + 4)*(AA - 4)*denom
    f0 = -alpha1*alpha2
    f1 = -f0 + alpha1 + alpha2 
    f2 = -alpha1 - alpha2 - 1
    E1 = EllipticCurve([0,f2,0,f1,f0])
    
    beta1 = (AA - 4*A + 4)*(BB - 4)*denom
    beta2 = (BB - 4*B + 4)*(AA - 4)*denom
    g0 = -beta1*beta2
    g1 = -g0 + beta1 + beta2 
    g2 = -beta1 - beta2 - 1
    E2 = EllipticCurve([0,g2,0,g1,g0])

    invariants = [type1_invariants,new_type1_invariants]
    matrices = []
    trafos = [special_trafo]
    codomain = [E1,E2]

    return [invariants, matrices, trafos, codomain]

def evaluate_splitting_map(P, type1_invariants):
    #compute the image of P on E1 and E2.
    # note: the second map is derived from map 1 by symmetry (A = -A, B = -B). 
    # this is not double-checked yet!
    [A,B,C,E] = type1_invariants
    assert C == 1,  "the invariants do not correspond to a splitting"
    [X0,X1,X2,X3] = P 

    AB = A*B

    P1x = E * (AB - A - A - B - B) * (X0 + X2)*(X0 + X2) + (4*E)*(X0 + X2)*X1 + 4*X3*(X0 - X1 + X2)
    P1z = E * (AB - A - A - B - B + 4) * (X0 - X2)*(X0 - X2)

    P2x = E * (AB + A + A + B + B) * (X0 + X2)*(X0 + X2) + (4*E)*(X0 + X2)*X1 + 4*X3*(X0 - X1 + X2)
    P2z = E * (AB + A + A + B + B + 4)* (X0 - X2)*(X0 - X2)

    return [[P1x,P1z],[P2x,P2z]]


def apply_splitting(P, phi):
    #ouput is a product of points on the Kummer lines
    [invariants, _, trafos, _] = phi
    [type1_invariants, new_type1_invariants] = invariants
    [special_trafo] = trafos
    P = apply_special_trafo(P, type1_invariants, special_trafo)
    P = evaluate_splitting_map(P, new_type1_invariants)
    return P
