from jacobian_isogenies.gluing_isogeny_jacobian import *
from jacobian_isogenies.formulae_jacobian import *
from jacobian_isogenies.doubling_jacobian import *
from sage.all import EllipticCurve


def splitting_jacobian_isogeny(ker4, phi):
    """
    input is a 4-torsion kernel lying above the next 2-isogeny, 
    and data of the previous (2,2)-isogeny
    We assume that the next (2,2)-kernel is of the form <((x-alpha),0), x*(x-beta, 0)>
    """
    invariants = phi[0]
    type2_invariants = invariants[-1]
    [Tp1,Tp2] = ker4
    ker2 = [double(Tp1,1,phi),double(Tp2,1,phi)]
    
    type1_invariants, trafo = find_special_trafo(type2_invariants,ker2,Tp1)
    [A,B,C,E] = type1_invariants
    assert C == 1, "the isogeny is not split"

    #splitting maps as in the paper (Prop. 4.2), we ignore twists since we cannot recover the y-coordinate anyways
    #W-points of the elliptic curves are alpha1,alpha2,1 
    # and beta1,beta2,1 respectively
    denom = 1/((A-2)*(B-2))
    alpha1 = (A+2)*(B-2)*denom
    alpha2 = (B+2)*(A-2)*denom
    f0 = -alpha1*alpha2
    f1 = -f0 + alpha1 + alpha2 
    f2 = -alpha1 - alpha2 - 1
    E1 = EllipticCurve([0,f2,0,f1,f0])
    
    denom = 1/(-f0)
    beta1 = alpha2*denom
    beta2 = alpha1*denom
    g0 = -beta1*beta2
    g1 = -g0 + beta1 + beta2 
    g2 = -beta1 - beta2 - 1
    E2 = EllipticCurve([0,g2,0,g1,g0])

    invariants = [type2_invariants,type1_invariants]
    codomain = [E1,E2]

    return [invariants, trafo, codomain]

def JacToKum(P, type1_invariants):
    [A,B,C,E] = type1_invariants
    assert C == 1, "this is not a split surface"
    [a0,a1,a2,b0,b1] = P
    y1y2 = b0*(b0*a2-b1*a1) + b1*b1*a0
    [X0,X1,X2] = [a2, -a1, a0]
    #note: we use C = 1 in the gluing step! general formula commented out
    #phi = E*(C*X0**2*X1 - 2*(A*C + B)*X0**2*X2 + (A*B + C + 1)*X0*X1*X2 - 2 * (A + B)*X0*X2**2 + X1*X2**2)
    phi = E*(X0*X0*(X1 - 2*(A + B)*X2) + (A*B + 2)*X0*X1*X2 + (-2*(A+B)*X0 + X1)*X2*X2)    
    X3 = (phi - 2*y1y2*X0*X0)/(X1*X1-4*X0*X2)
    #assert (C**2*E**2)*X0**4 + (-2*A*B*C*E**2 - 2*C**2*E**2 - 2*C*E**2)*X0**3*X2 + (4*A*C*E**2 + 4*B*C*E**2)*X0**2*X1*X2 + (-4*C*E**2)*X0*X1**2*X2 + (A**2*B**2*E**2 - 4*A**2*C*E**2 - 2*A*B*C*E**2- 2*A*B*E**2 - 4*B**2*E**2 + C**2*E**2 + 4*C*E**2 + E**2)*X0**2*X2**2 + (4*A*C*E**2 + 4*B*E**2)*X0*X1*X2**2 + (-2*A*B*E**2 - 2*C*E**2 - 2*E**2)*X0*X2**3 + (E**2)*X2**4 + (-2*C*E)*X0**2*X1*X3 + (4*A*C*E + 4*B*E)*X0**2*X2*X3 + (-2*A*B*E - 2*C*E - 2*E)*X0*X1*X2*X3 + (4*A*E + 4*B*E)*X0*X2**2*X3 + (-2*E)*X1*X2**2*X3 + X1**2*X3**2 - 4*X0*X2*X3**2 == 0, P
    return [X0,X1,X2,X3]



def evaluate_splitting_map(P, type1_invariants):
    #compute the image of P on E1 and E2.
    # note: the second map is derived from map 1 by symmetry (A = -A, B = -B). 
    # this is not double-checked yet!
    [A,B,C,E] = type1_invariants
    assert C == 1,  "the invariants do not correspond to a splitting"
    [X0,X1,X2,X3] = JacToKum(P, type1_invariants)

    P1x = E * (A*B - A - A - B - B) * (X0 + X2)*(X0 + X2 + 4*E*X1) + 4*X3*(X0 - X1 + X2)
    P1z = E * (A - 2)*(B - 2)* (X0 - X2)*(X0 - X2)

    P2x = E * (A*B + A + A + B + B) * (X0 + X2)*(X0 + X2 + 4*E*X1) + 4*X3*(X0 - X1 + X2)
    P2z = E * (A + 2)*(B + 2)* (X0 - X2)*(X0 - X2)

    return [[P1x,P1z],[P2x,P2z]]


def apply_splitting(P, phi):
    #ouput is a product of points on the Kummer lines
    [invariants, trafo, _] = phi
    [type2_invariants, type1_invariants] = invariants
    P = apply_Mumford_transformation(P, trafo)
    P = evaluate_splitting_map(P, type1_invariants)
    return P