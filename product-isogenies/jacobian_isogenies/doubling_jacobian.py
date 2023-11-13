
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


def double_type2(P, mult, phi):
	"""
	compute 2**mult*P, where P is a point on the Jacobian of 
	y**2 = (x**2-1)(x**2-A)(E*x**2-B*x+C)
	and [A,B,C,E] are the type2_invariants
	obtained from the previous isogeny phi
	"""

	[A,B,C,E] = phi[0][-1]
	[t1,AC] = phi[2][-1]
	# Precomputations
	#t1 = (A+1)*E - C
	#AC = A*C
	for i in range(mult):
		[a0,a1,a2,b0,b1] = P
		if a2 == 0:
			print("not implemented");
		# Precomputations
		a1b1 = a1*b1;
		a02 = a0*a0;
		a12 = a1*a1;
		b02 = b0*b0;
		a0b1 = a0*b1;
		a1b0 = a1*b0;
		# Composition step: 2*P = (a**2, bb), where bb = b + a * (e1x+e0)
		Edeni = 1/((a0+a0)*(b1*(a1b0 - a0b1) - b02));
		E0 = a02*(E*(a12*(a1b1 + b0+b0+b0) - (a0+a0)*(a1b1 + a1b1 + b0))+ B*((a12-A-1)*b1  + a1b0 + a1b0 - a0b1 - a0b1) -t1*(a1b1+b0)) + (AC-b02)*(a1b1-b0);
		E1 = a02*((E+E)*(-a1b0-a1b0-a1b0  + a0b1) - ((E+E+E)*a1 +B+B)*a1b1)  + a0*((E+E+E+E)*a12*a1b0 +t1*(a0b1-a1b0-a1b0) - (a0+a0+A+1-a12-a12-a12)*B*b0) + b1*(AC- b02);
		e0 = E0*Edeni;
		e1 = E1*Edeni;
		b00 = a0*e0+b0;
		b11 = a1*e0+a0*e1+b1;
		b22 = a1*e1+e0;
		b33 = e1;

		den = b33*b33 - E;
		if den == 0:
			a0p = (-AC + b00*b00)/(a02*(B + 2*b22*b33));
			a1p = 1;
			a2p = 0;
			b0p = (b33*a0p - b22)*a0p*a0p + b11*a0p - b00;
			b1p = 0;
		else:
			denom = 1/(a02*den)
			a0p = (b00*b00-AC)*denom;
			a1p = -a1-a1 + (B  + (b22+b22)*b33)*denom*a02;
			a2p = 1;
			b1p = -b33*(a1p*a1p - a0p) + b22*a1p - b11;
			b0p = a0p*(-b33*a1p + b22) - b00;
		P = [a0p,a1p,a2p,b0p,b1p]
	return P


def double(P,mult,phi=None):
	if len(P) == 2:
		return double_prod(P,mult)
	else:
		return double_type2(P,mult,phi)

