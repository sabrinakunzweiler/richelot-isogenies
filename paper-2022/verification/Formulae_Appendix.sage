"""
Code accompanying the paper "Efficient Computation of (2^n,2^n)-Isogenies"
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

This code can be used to verify the formulas deduced in Appendix A of the paper.
"""


#Type 1 Equation:
K.<u,A,B,C,E,a0,a1,b0,b1> = QQ[]
S.<x> = K[]
F = E*x*(x^2-A*x+1)*(x^2-B*x+C)
Ap = C
Bp = 2/E
Cp = (B-A*C)/(E*(1-C))
Ep = (A-B)/(E*(1-C))
Fnew =  (x^2-1)*(x^2-Ap)*(Ep*x^2-Bp*x+Cp)

#Richelot correspondence:
rel1 =  (u^2-B*u+1)*x^2 + 2*(C-1)*u*x - C*u^2+B*u-C
rel2 = (A-B)*u*x^3 - ((A-B)*u^2+2*(1-C)*u)* x^2
+ (2*(1-C)*u^2 - (A*C-B)*u)*x + (A*C-B)*u^2

#relations among Mumford coefficients
a = x^2+a1*x+a0
b = b1*x+b0
[q,r] = (F-b^2).quo_rem(a) #r must be zero
relations = [a(u)] + r.coefficients()
I = K.ideal(relations)
v = b0 + b1*u
s = -a1 - u
t = b0 + b1*s

print("Lemma A.1:", -b1*(a1*b0-a0*b1)+b0^2 == a.resultant(b))

print("Proposition A.4")
print("v1:", rel2(u=-a0, x = B/2) == 1/8*(4*C-B^2)*(A-B)*a0*(B+2*a0))
print("v2:", (rel2(u=-a0, x=1/x)*x^3)(0) == (B-A)*a0)

print("Lemma A.5:", a0*B^2 + (a0+1)*a1*B+(a0-1)^2+a1^2 == a.resultant(x^2-B*x+1))
#general check that Q_1=(s1,t1), Q2 =(inf_sign) are indeed on the curve
s = -a1-u
t = b1*s+b0
J = I + K.ideal(s^2-B*s+1)
s1 = B/2
t1 = (4*C-B^2)*(B-A)*s*(B-2*s)/(8*t)
e = (A-B)*s/t
print("General checks")
print("Q1 on curve:", (t1^2-Fnew(s1)).numerator().reduce(J) == 0)
print("square-root E", (Ep-e^2).numerator().reduce(J) == 0)

print("Proposition A.7:")
J = I + K.ideal([a1+B,a0-1])
u1 = B/2
v1 = (4*C-B^2)*(B-A)*u*(B-2*u)/(8*v)
print("case b0=0:", (v1/t1).numerator().reduce(J + K.ideal(u*t-s*v))
== -1*(v1/t1).denominator().reduce(J + K.ideal(u*t-s*v)) )
print("case b0!=0:",
(v1 - (4*C-B^2)*(B-A)/(4*b0)).numerator().reduce(J+K.ideal(u*t+s*v)) == 0)

#Section A.3
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

print("Lemma A.8:", (aQ.resultant(aP).numerator() + (C-1)^2*(u-s)^2*(4*C-B^2)*(1-u*s)^2).reduce(I) == 0)

print("Proposition A.10:")
J = I + K.ideal([a0-1])
rel1_s = (B+a1)*x^2 - 2*(C-1)*x - (B+a1*C)
print(all([c.reduce(J) == 0 for c in (rel1 +u*rel1_s).coefficients()]))


d1 = (B*b0 - 2*(a1*b0-b1))*(a1+B) -4*b0*(C-1)
d0 = (B*(b1-a1*b0)+ 2*C*b0)*(a1+B) - 2*B*b0*(C-1)
#nonzero terms
nz = [B+a1,4*C-B^2, C-1, A+a1, -a1*b0*b1+b0^2+b1^2]
print("bP1-bQ1:", all([(bP1-bQ1).numerator().reduce(J) == -nz[0]^2*nz[2]*nz[3]*(2*u+a1)*d1, (bP1-bQ1).denominator().reduce(J) == (v*t*nz[0]^4).reduce(J)]))
print("bP0-bQ0:", all([(bP0-bQ0).numerator().reduce(J) == nz[0]^2*nz[2]*nz[3]*(2*u+a1)*d0, (bP1-bQ1).denominator().reduce(J) == (v*t*nz[0]^4).reduce(J)]))
xhat = (bQ0-bP0)/(bP1-bQ1)
print("xhat:", (xhat - d0/d1).numerator().reduce(J) == 0)
yhat = bP1*xhat + bP0
print("yhat:", (yhat - nz[1]*nz[2]*nz[3]/d1).numerator().reduce(J) == 0)
print("xhat is a root of aP:", rel1(xhat).numerator().reduce(J) == 0)

print("check that d1 nonzero (by contradiction):")
print("if d1=d0=0, ")
J1 = J + K.ideal([d0,d1])
print("then b1=b0=0 (contradiction):", all([prod(nz)^2*b0 in J1, prod(nz)^2*b1 in J1]))
print("if d0 nonzero, then bP1=bQ1=0.", True) #geometric argument
J2 = J + K.ideal([d1, (bP1+bQ1).numerator().reduce(J)])
print("then 0=1 (contradiction):", prod(nz)^2 in J2)
