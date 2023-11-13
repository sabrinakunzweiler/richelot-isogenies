from jacobian_isogenies.formulae_jacobian import *
from jacobian_isogenies.doubling_jacobian import *


def jacobian_isogeny(ker4, phi):
	"""
	ker4 are 4-torsion points lying above the (2,2)-isogeny that is computed
	phi is the previous isogeny 
	"""
	invariants = phi[0]
	type2_invariants = invariants[-1]
	ker2 = [double(ker4[0],1,phi), double(ker4[1],1,phi)]
	type1_invariants, special_trafo = find_special_trafo(type2_invariants, ker2, ker4[0])
	codomain_type2_invariants = type1_isogeny(type1_invariants)
	invariants = [type2_invariants, type1_invariants,codomain_type2_invariants]
	aux_iso = auxiliary_products(type1_invariants)
	aux_doub = auxiliary_doubling(codomain_type2_invariants)
	aux = [aux_iso,aux_doub]

	return [invariants, special_trafo, aux]


def apply_jacobian_isogeny(P, phi):
	"""
	Push point P through the isogeny phi. 
	This means applying applying a transformation to Type 1
	and then applying the isogeny to Type 2
	"""
	[invariants, trafo, aux] = phi
	type1_invariants = invariants[1]
	aux_iso = aux[0]
	P = apply_Mumford_transformation(P,trafo)  
	P = apply_type1_isogeny(P,type1_invariants,aux_iso)

	return P