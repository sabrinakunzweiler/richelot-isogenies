from kummer_isogenies.doubling_kummer import *

def find_elliptic_transformation(kernel2,kernel4):
    """
    we assume that kernel4 = [[P1,P2],[Q1,Q2]] where P1,Q1 are 4-torsion points on E1 and P2,Q2 are 4-torsion points on E2
    similarly for kernel2
    Output:
    transformations [a_r,b_r,e_r] and [a_s,b_s,e_s] for the curves E1 and E2 so that the image curves 
    have the correct form for our gluing formulae.
    """
    
    E1 = kernel2[0][0].curve()
    E2 = kernel2[0][1].curve()

     #compute r1,r2,r3 Weierstrass points of E1, and s1,s2,s3 Weierstrass points of E2.
    r3 = kernel2[0][0][0]
    r2 = kernel2[1][0][0]
    r1 = - E1.a2() - r3 - r2 #we assume standard Weierstrass model with a1 = a3 = 0
    #assert r1*r2*r3 == 0, 'something probably went wrong'
    s3 = kernel2[0][1][0]
    s2 = kernel2[1][1][0]
    s1 = - E2.a2() - s3 - s2
    
    #avoid inversions
    large_denom = 1/(2*kernel4[0][0][1]*kernel4[0][1][1]*(r3 - r1)*(r2 - r1)*(s3 - s2)*(s3 - s1)*(s2 - s1)*(r3 - r2))

    #extract canonical square-roots for the transformations
    rts1_prod = kernel4[0][0][0] - r3
    rts1_sum_inv = rts1_prod*kernel4[0][1][1]*(r3 + r3 - r1 - r1)*(r2 - r1)*(s3 - s2)*(s3 - s1)*(s2 - s1)*(r3 - r2)*large_denom
    sqr3r2 = (rts1_prod + r3 - r2)*rts1_sum_inv
    #assert sqr3r2**2 == r3-r2, "sqr3r2"
    
    rts2_prod = kernel4[0][1][0] - s3
    rts2_sum_inv = rts2_prod*kernel4[0][0][1]*(r3 +r3 - r1 - r1)*(r2 - r1)*(s3 - s2)*(s3 - s1)*(s2 - s1)*(r3 - r2)*large_denom
    sqs3s2 = (rts2_prod + s3 - s2)*rts2_sum_inv
    #assert sqs3s2**2 == s3-s2, "sqs3s2"
    
    #transformation a_r*x+b_r, e_r*y on E1
    term1 = (r2 - r3)*s1 + (r3 - r1)*s2  + (r1 - r2)*s3
    den1 = kernel4[0][0][1]*kernel4[0][1][1]*(s3 - s1)*(s2 - s1)*(r3 - r2)*large_denom
    a_r = (term1+term1)*den1
    b_r = 1 - r1*a_r
    e_r = sqs3s2*den1
    #transformation a_s*x+b_s, e_s*y on E2
    den2 = kernel4[0][0][1]*kernel4[0][1][1]*(r3 - r1)*(r2 - r1)*(s3 - s2)*large_denom
    a_s = -(term1+term1)*den2
    b_s = 1 - s1*a_s
    e_s = sqr3r2*den2

    #new Weierstrass points are 1,alpha1,alpha2 and 1,beta1,beta2
    alpha1 = a_r*r2 + b_r 
    #alpha2 = a_s*s2 + b_s 
    beta1 = a_r*r3 + b_r 
    #beta2 = a_s*s3 + b_s 
    #assert a_r*r1 + b_r == 1
    #assert a_s*s1 + b_s == 1
    #parameters A and B
    denom = 1/((alpha1 - 1)*(beta1 - 1)*term1)
    A = 2*(alpha1 + 1)*(beta1 - 1)*term1*denom
    B = 2*(beta1 + 1)*(alpha1 - 1)*term1*denom
    E = (alpha1 - 1)*(beta1 - 1)*denom
    #assert alpha1 == (A+2)/(A-2), "alpha1"
    #assert alpha2 == (A-2)/(A+2), "alpha2"
    #assert beta1 == (B+2)/(B-2), "beta1"
    #assert beta2 == (B-2)/(B+2), "beta2"
    #assert E == 8**2*e_r**2/((A-2)*(B-2)*a_r**3), "E1"
    #assert E == - 8**2*e_s**2/((A+2)*(B+2)*a_s**3), "E2"
    
    return [[A,B,1,E], [[a_r,b_r,e_r],[a_s,b_s,e_s]]]
    
def elliptic_transformation(P, trafo):
        #input: trafo = [a,b,e] defining transformation of the form x -> a*x + b and y -> e*y on the level of points!
        [a,b,e] = trafo
        if P != 0:
            return [a*P[0] + b, e*P[1]]
        else:
            return 0

def elliptic_transformation_couple(couple,trafos):
    #input is a pair of points on E1 x E2, we apply different transformations
    image_1 = elliptic_transformation(couple[0],trafos[0])
    image_2 = elliptic_transformation(couple[1],trafos[1])

    return [image_1,image_2]

def type1_add(type1_invariants, D1,D2):
    #Input: *type1_invariants = [A,B,C,E] defining hyperelliptic curve CC: y**2 = E*x*(x**2-A*x+1)*(x**2-B*x+C)
    #       * D1, D2 are divisors given in Mumford presentation on Jac(C)
    #Output: * Mumford coordinates of D1 + D2 
    #Note: ad-hoc solution
    [a,b] = D1
    [c,d] = D2
    [A,B,C,E] = type1_invariants
    K = A.parent()
    R = K['x']; (x,) = R._first_ngens(1)
    F = E*x*(x*x-A*x+1)*(x*x-B*x+C)
    #assert (F-b**2) % a == 0, "a,b incorrect"
    #assert (F-d**2) % c == 0, "c,d incorrect"
    a_comp = a*c
    den = (a[1]*c[0] - a[0]*c[1])*(a[1]-c[1]) + (a[0]-c[0])**2
    coef1 = (a[1]-c[1])*(b[0] - d[0]) + (a[0]-c[0])*(d[1]-b[1])
    coef0 = (a[1]-c[1])*(c[1]*(b[0] - d[0]) - c[0]*(b[1] - d[1])) + (a[0]-c[0])*(d[0] - b[0])
    b_comp = b + a*(coef1*x+coef0)/den
    [new_a,rem] = (F-b_comp*b_comp).quo_rem(a_comp)
    new_b = -b_comp % new_a
    return [new_a,new_b]
    
def JacToKum(P, type1_invariants):
    [A,B,C,E] = type1_invariants
    assert C == 1, "this is not a split surface"
    a = P[0]
    b = P[1]
    y1y2 = b[0]*b[0]*a[2]-b[0]*b[1]*a[1] + b[1]*b[1]*a[0]
    [X0,X1,X2] = [a[2], -a[1], a[0]]
    #note: we use C = 1 in the gluing step! general formula commented out
    #phi = E*(C*X0**2*X1 - 2*(A*C + B)*X0**2*X2 + (A*B + C + 1)*X0*X1*X2 - 2 * (A + B)*X0*X2**2 + X1*X2**2)
    phi = E*(X0*X0*(X1 - 2*(A + B)*X2) + (A*B + 2)*X0*X1*X2 + (-2*(A+B)*X0 + X1)*X2*X2)    
    X3 = (phi - 2*y1y2*X0*X0)/(X1*X1-4*X0*X2)
    #assert (C**2*E**2)*X0**4 + (-2*A*B*C*E**2 - 2*C**2*E**2 - 2*C*E**2)*X0**3*X2 + (4*A*C*E**2 + 4*B*C*E**2)*X0**2*X1*X2 + (-4*C*E**2)*X0*X1**2*X2 + (A**2*B**2*E**2 - 4*A**2*C*E**2 - 2*A*B*C*E**2- 2*A*B*E**2 - 4*B**2*E**2 + C**2*E**2 + 4*C*E**2 + E**2)*X0**2*X2**2 + (4*A*C*E**2 + 4*B*E**2)*X0*X1*X2**2 + (-2*A*B*E**2 - 2*C*E**2 - 2*E**2)*X0*X2**3 + (E**2)*X2**4 + (-2*C*E)*X0**2*X1*X3 + (4*A*C*E + 4*B*E)*X0**2*X2*X3 + (-2*A*B*E - 2*C*E - 2*E)*X0*X1*X2*X3 + (4*A*E + 4*B*E)*X0*X2**2*X3 + (-2*E)*X1*X2**2*X3 + X1**2*X3**2 - 4*X0*X2*X3**2 == 0, P
    return [X0,X1,X2,X3]


def pullback1(P):
    if P!=0:
        [u1,v1] = P
        K = u1.parent()
        R = K['x']; (x,) = R._first_ngens(1)
        u11_inv = 1/(u1-1)
        return [x**2 -2*(u1+1)*u11_inv*x + 1, 4*v1*u11_inv*((1 + 4*u11_inv)*x - 1)]
    else:
        return 0

def pullback2(Q):
    if Q!=0:
        [u2,v2] = Q
        K = u2.parent()
        R = K['x']; (x,) = R._first_ngens(1)
        u21_inv = 1/(u2-1)
        return [x**2 + 2*(u2+1)*u21_inv*x + 1, 4*v2*u21_inv*((1 + 4*u21_inv)*x  + 1)]
    else:
        return 0

def apply_gluing(P,gluing):
    [invariants, trafos] = gluing
    type1_invariants = invariants[-1]
    P = elliptic_transformation_couple(P,trafos)
   
    image_1 = pullback1(P[0])
    image_2 = pullback2(P[1])
    if image_1 == 0:
        image = image_2
    elif image_2 == 0:
        image = image_1 
    else:
        image = type1_add(type1_invariants, image_1,image_2)
    image = JacToKum(image,type1_invariants)
    return image

    
def gluing_kummer_isogeny(kernel4):
    #Input: kernel = [[P1,P2],[Q1,Q2]] generates a 4-maximal isotropic
    #       on a product of elliptic curves E1 x E2
    #       both E1 and E2 are assumed to be in Montgomery Form: y**2 = x**3 + Ai*x**2 +x
    #Output: *The type-1 invariants [A,B,C,E] defining a hyperelliptic curve C: y**2 = E*x*(x**2-A*x+1)*(x**2-B*x+C)
    #        so that there is a is (2,2)-isogeny phi: E1xE2 -> Jac(C) with ker(phi) = 2**(n-1)*kernel
    #        * the transformations to the correct form of elliptic curve to apply the gluing map
    
    [P1, P2] = kernel4[0]
    [Q1, Q2] = kernel4[1]

    kernel2 = [[2*P1,2*P2],[2*Q1,2*Q2]]

    [type1_invariants, trafos] = find_elliptic_transformation(kernel2,kernel4)
    
    return [[type1_invariants], trafos]



