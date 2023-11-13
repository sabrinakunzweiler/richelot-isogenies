from kummer_isogenies.formulae_kummer import *

def double_montgomery(P, mult):
	"""
	compute 2**mult*P, where P is a point on an elliptic curve in Montgomery form.
	"""
	E = P.curve()
	A = E.a2()
	[x1,y1,z1] = P
	for i in range(mult):
		m = ((x1+x1+x1)*x1 + A*(x1+x1) + 1)/(y1 + y1)
		m2 = m*m
		x2 = m2 - (A+x1+x1)
		y1 = ((x1+x1+x1+A) - m2)*m - y1
		x1 = x2
	return E([x1, y1,1])

def double_prod(P, mult):
	"""
	compute 2**mult*P, where P = (P1,P2) on elliptic product
	"""
	[P1,P2] = P
	return [double_montgomery(P1,mult),double_montgomery(P2,mult)]


def double_kummer(P,mult,phi):
	"""
	compute 2**mult*P, where P is a point a type-1 Kummer surface
	phi contains the data of the previous (2,2)-isogeny (we reuse the matrix M)
	"""
	matrices = phi[1]
	[M,L] = matrices
	for i in range(mult):
		P = apply_Richelot_type1(P,L)
		P = apply_Richelot_type2(P,M)
	return P


def double(P,mult,phi=None):
	if len(P) == 2:
		return double_prod(P,mult)
	else:
		assert phi, "need data of the previous isogeny"
		return double_kummer(P,mult,phi)

