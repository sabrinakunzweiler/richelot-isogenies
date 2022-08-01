"""
Code accompanying the paper "Efficient Computation of (2^n,2^n)-Isogenies"
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

This code can be used to verify the formula in Corollary 2.8.
"""


R.<r1,r2,r3,r4> = PolynomialRing(QQ)
s1 = r1 + r2 + r3 + r4
s2 = r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4
s3 = r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4
s4 = r1*r2*r3*r4
b1 = -r1^2
b2 = -r2^2

print(r1*r2 == (s1*s3*b1*b2 + (s4-b1*b2)^2) / (b1*b2*s1^2 + (s4-b1*b2)*(s2-b1-b2)))
