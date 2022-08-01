"""
Code accompanying the paper "Efficient Computation of (2^n,2^n)-Isogenies"
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

This code can be used to verify the formulas deduced in Section 4 of the paper.
"""

def Richelot(G, delta):
    #input: list of 3 quadratic polynomials defining the isogeny, and the value of delta (neq 0)
    #output: list of three quadratic polynomials [h1,h2,h3] s.t. y^2=h1*h2*h3 is the new equation
    assert len(G) == 3
    Gd = [g.derivative() for g in G]
    H = [(Gd[(i+1)%3]*G[(i+2)%3]-Gd[(i+2)%3]*G[(i+1)%3])/delta for i in range(3)]
    return H

#Type 1 Equation:
R.<A,B,C,E,u,v> = QQ[]
S.<x> = Frac(R)[]
F = E*x*(x^2-A*x+1)*(x^2-B*x+C)
# special kernel of the isogeny is G
G = [E*x,x^2-A*x+1, (x^2-B*x+C)]


print("Proposition 4.2 ")
delta = -E*(1-C)
H = Richelot(G,delta);
Ap = C
Bp = 2/E
Cp = (B-A*C)/(E*(1-C))
Ep = (A-B)/(E*(1-C))
print(prod(H)*(1-C)^2 == (x^2-1)*(x^2-Ap)*(Ep*x^2-Bp*x+Cp))

print("Proposition 4.3")
P.<up> = Frac(R)[]
rel1 = (G[0](u)*H[0](up)+G[1](u)*H[1](up))*(1-C);  # note that here multiplication by a constant is ok
rel2 = (G[0](u)*H[0](up))*(u-up)*(1-C);  #rel2=(1-C)*v'*v
print(rel1 == -(u^2-B*u+1)*up^2 - 2*(C-1)*u*up + C*u^2-B*u+C)
print(rel2 == (A-B)*u*up^3 - ((A-B)*u^2+2*(1-C)*u)* up^2 + (2*(1-C)*u^2 - (A*C-B)*u)*up + (A*C-B)*u^2 )

print("Lemma 4.5")
aP1 = 2*(C-1)*u/(u^2-B*u+1)
aP0 = (-C*u^2+B*u-C)/(u^2-B*u+1)
print(rel1/(-u^2+B*u-1) == up^2 + aP1*up + aP0)
bP1 = u*(1-C)*(u^2-A*u+1)/(u^2-B*u+1)^2 * (2*u^3-B*u^2 + (-B^2+4*C-2)*u+B)
bP0 = -u*(1-C)*(u^2-A*u+1)/(u^2-B*u+1)^2 * (B*u^3+(-B^2+2*C)*u^2 - B*u+2*C)
print(rel2 % rel1 == bP1*up + bP0)


print("Theorem 4.7")
K.<A,B,C,E,u,a0,a1,b0,b1> = QQ[]
R.<x> = K[]
#Relations among the elements
# 1) u is a root of a(x) = x^2+a1*x+a0
# 2) a0,a1,b0,b1 describe a divisor on the curve y^2 = x(x^2-Ax+1)(x^2-Bx+C)
rela = u^2 + a1*u + a0
F = E*x*(x^2-A*x+1)*(x^2-B*x+C)
b = b1*x+b0
a = x^2+a1*x+a0
[q,r] = (F-b^2).quo_rem(a) #r must be zero
relations = [rela] + r.coefficients()
I = K.ideal(relations)
v = b0 + b1*u
s = -a1 - u
t = b0 + b1*s
#expressions for aP, bP already verified above
aP1 = 2*(C-1)*u/(u^2-B*u+1)
aP0 = (-C*u^2+B*u-C)/(u^2-B*u+1)
bP1 = u*(1-C)*(u^2-A*u+1)/(u^2-B*u+1)^2/v * (2*u^3-B*u^2 + (-B^2+4*C-2)*u+B)
bP0 = -u*(1-C)*(u^2-A*u+1)/(u^2-B*u+1)^2/v * (B*u^3+(-B^2+2*C)*u^2 - B*u+2*C)
aP = x^2+aP1*x+aP0
bP = bP1*x +bP0

aQ1 = 2*(C-1)*s/(s^2-B*s+1)
aQ0 = (-C*s^2+B*s-C)/(s^2-B*s+1)
bQ1 = +s*(1-C)*(s^2-A*s+1)/(s^2-B*s+1)^2/t * (2*s^3-B*s^2 + (-B^2+4*C-2)*s+B)
bQ0 = -s*(1-C)*(s^2-A*s+1)/(s^2-B*s+1)^2/t * (B*s^3+(-B^2+2*C)*s^2 - B*s+2*C)
aQ = x^2+aQ1*x+aQ0
bQ = bQ1*x +bQ0

a00 =  a0*B^2 + (a0*a1 + a1)*B*C + (a0^2 + a1^2 - 2*a0 + 1)*C^2
a11 =  4*a0*B*C + (2*a0*a1 + 2*a1)*C^2 + (-4*a0)*B + (-2*a0*a1 - 2*a1)*C
a22 =  (-2*a0)*B^2 + (-a0*a1 - a1)*B*C + 4*a0*C^2 + (-a0*a1 - a1)*B + (-2*a0^2 - 2*a1^2 - 4*a0 - 2)*C + 4*a0
a33 =  (-4*a0)*B*C + 4*a0*B + (-2*a0*a1 - 2*a1)*C + 2*a0*a1 + 2*a1
aden =  a0*B^2 + (a0*a1 + a1)*B + a0^2 + a1^2 - 2*a0 + 1
ap = (a00 + a11*x + a22*x^2 + a33*x^3+aden*x^4)/aden
#need to use relations to find that aP*aQ = a. We multiply by aden to have integral coefficients
acomp = (aP*aQ).coefficients()
print("representation for a':")
print("a0':", K(a00-acomp[0].numerator()).reduce(I) == 0)
print("a1':", K(a11-acomp[1].numerator()).reduce(I) == 0)
print("a2':", K(a22-acomp[2].numerator()).reduce(I) == 0)
print("a3':", K(a33-acomp[3].numerator()).reduce(I) == 0)
print("a4':", all([c.denominator().reduce(I) == aden for c in acomp[:3]]))

b00 =  (a0*a1*b0 - a0^2*b1)*A*B + (a0^2*b0 + a1^2*b0 - a0*a1*b1 - a0*b0)*A*C + (a0*a1^2*b0 - a0^2*a1*b1 - a0^2*b0 + a0*b0)*B + (a0^2*a1*b0 + a1^3*b0 - a0^3*b1 - a0*a1^2*b1 - 2*a0*a1*b0 + 2*a0^2*b1 + a1*b0 - a0*b1)*C
b11 =  a0*b0*A*B + (2*a0*a1*b0 - a0^2*b1 + a1*b0 - a0*b1)*A*C + (-2*a0*a1*b0 + 2*a0^2*b1)*A + (a0*a1*b0 - a0^2*b1 + a0*b1)*B + (2*a0*a1^2*b0 - 2*a0^2*a1*b1 - a0^2*b0 + a1^2*b0 + b0)*C - 2*a0*a1^2*b0 + 2*a0^2*a1*b1 + 2*a0^2*b0 - 2*a0*b0
b22 =  (-a0*a1*b0 + a0^2*b1)*A*B + 2*a0*b0*A*C + (-a0^2*b0 - a1^2*b0 + a0*a1*b1 - a0*b0)*A + (-a0*a1^2*b0 + a0^2*a1*b1 + a0^2*b0 - a0*b0)*B + (2*a0*a1*b0 - 2*a0^2*b1 + 2*a0*b1)*C - a0^2*a1*b0 - a1^3*b0 + a0^3*b1 + a0*a1^2*b1 - a1*b0 - a0*b1
b33 =  (-a0*b0)*A*B + (-a0^2*b1 - a1*b0 + a0*b1)*A + (-a0*a1*b0 + a0^2*b1 - a0*b1)*B - a0^2*b0 - a1^2*b0 + 2*a0*b0 - b0
bden =  -1*(a0 - 1) * (-a1*b0*b1 + a0*b1^2 + b0^2)
bp = (b33*x^3+b22*x^2+b11*x+b00)/bden
print("representation for b':")
print("b' = bP (mod aP): ", all([c.numerator().reduce(I) == 0 for c in ((bp-bP)%aP).coefficients()]))
print("b' = bQ (mod aQ): ", all([c.numerator().reduce(I) == 0 for c in ((bp-bQ)%aQ).coefficients()]))

Ap = C
Bp = 2/E
Cp = (B-A*C)/(E*(1-C))
Ep = (A-B)/(E*(1-C))
Fp = (x^2-1)*(x^2-Ap)*(Ep*x^2-Bp*x+Cp)
print("b'^2 = f (mod a'): ", all([c.numerator().reduce(I) == 0 for c in ((Fp-bp^2)%ap).coefficients()]))
