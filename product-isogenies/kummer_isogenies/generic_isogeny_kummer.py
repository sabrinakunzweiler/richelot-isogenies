from kummer_isogenies.formulae_kummer import *
from kummer_isogenies.doubling_kummer import *


def kummer_isogeny(ker4,phi):
    """
    ker4 are 4-torsion points lying above the (2,2)-isogeny that is computed
    phi is the previous isogeny 
    """

    invariants = phi[0]
    type1_invariants = invariants[-1]
    [Tp1,Tp2] = ker4
    new_Tp1 = double(Tp1,1,phi)
    new_Tp2 = double(Tp2,1,phi)
    ker2 = [new_Tp1,new_Tp2]
    [new_type1_invariants, special_trafo] = find_special_trafo(type1_invariants, ker2, Tp1)
    trafo_12 = compute_trafo_12(new_type1_invariants) #this is for transformation from type 1 to 2
    new_type2_invariants = invariant_type1_type2(new_type1_invariants)
    M = Richelot_type2_matrix(new_type2_invariants)
    #data below is needed for the doubling
    codomain_type1_invariants = Richelot_type2_invariants(new_type2_invariants)
    L = Richelot_type1_matrix(codomain_type1_invariants) #this is the dual to M
    invariants = [type1_invariants, new_type2_invariants,codomain_type1_invariants]
    matrices = [M,L]
    trafos = [special_trafo,trafo_12]

    return [invariants, matrices, trafos]

def apply_kummer_isogeny(P, phi):
    [invariants, matrices, trafos] = phi
    M = matrices[0]
    [type1_invariants, new_type2_invariants,_] = invariants
    [special_trafo,trafo_12] = trafos
    P = apply_special_trafo(P, type1_invariants, special_trafo)
    P = apply_trafo_type1_type2(P, trafo_12)
    P = apply_Richelot_type2(P, M)
    return P

def first_kummer_isogeny(ker, phi, ker2):
    """
    Here we also require the 2-subgroup
    """
    invariants = phi[0]
    type1_invariants = invariants[-1]
    [Tp1,Tp2] = ker
    [new_type1_invariants, special_trafo] = find_special_trafo(type1_invariants, ker2, Tp1)
    new_type2_invariants = invariant_type1_type2(new_type1_invariants)
    trafo_12 = compute_trafo_12(new_type1_invariants) #this is for transformation from type 1 to 2
    M = Richelot_type2_matrix(new_type2_invariants)
    #data below is needed for the doubling
    codomain_type1_invariants = Richelot_type2_invariants(new_type2_invariants)
    L = Richelot_type1_matrix(codomain_type1_invariants) #this is the dual to M
    invariants = [type1_invariants, new_type2_invariants,codomain_type1_invariants]
    matrices = [M,L]
    trafos = [special_trafo, trafo_12]

    return [invariants, matrices, trafos]
