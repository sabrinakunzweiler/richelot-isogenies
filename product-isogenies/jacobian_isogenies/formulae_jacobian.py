def find_special_trafo(type2_invariants,ker2,P4):
    """
    INPUT: type2_invariants =[A,B,C,E] defining the hyperelliptic curve
    y**2 = (x**2-1)(x**2-A)(Ex**2-Bx+C)
    the 4-torsion kernel ker 4
    OUTPUT: [alpha,beta,gamma,delta,epsilon] defining the coordinate transformation
    from Proposition 4.1 to a Type-1 equation
    [Ap,Bp,Cp,Ep] coeffcients of the Type-1 equation
    """

    [A,B,C,E] = type2_invariants;
    [g10,g11,g12] = ker2[0][:3]
    [g20,g21,g22] = ker2[1][:3]

    # proceed as in Proof of Proposition 4.5
    assert g10 == -g11-1  #check that indeed 1 is a root of g1
    alpha1 = 1
    alpha2 = g10

    """
    either -alpha1 or -alpha2 are a root of g2, which we call beta1,
    the other root is called gamma1
    : beta2 and gamma2 are the roots of (Ex**2-Bx+C)
    Note: We directly apply the transformation t1: x -> (x-alpha2)/(x-alpha1)
    """
    lc = (alpha1-alpha2)*(alpha1-alpha2);
    den_aux = 1/((g21+g21)*g11)
    assert g20 == g21 - alpha1
    beta1 = -g11*g21*g11*den_aux;
    beta2 = (g21-alpha1+alpha2)*g11*(den_aux+den_aux);
    gamma1 = -(alpha2+alpha2)*g21*(den_aux+den_aux);
    gamma2den = B+E*(g21-2);
    gamma2 = (B+E*(g21+g11))/gamma2den;
    lc = (lc+lc)*g21*g11*gamma2den;
    
    # Transformation (x-alpha2)/(x-alpha1) is applied  to the 4-torsion point P
    [a0,a1,a2,b0,b1] = P4;
    mu = a1*b0-a0*b1;
    deni = 1/((a1 + 1) + a0);
    a0p =  (alpha2*(a1 + alpha2) + a0)*deni;
    a1p =  (a1*g11 - (g10 + a0+ g10 + a0))*deni;
    aux1 = (deni*(alpha2 - 1))**2;
    b0p =  aux1 * (((b1*alpha2+b0) + (b0+b0)*alpha2) + mu*(2+alpha2+a1) - a0*b0);
    b1p =  aux1 * (a0*b0 - (b1 + b0+b0+b0) - (3 + a1)*mu);

    # compute sqrt(beta1*beta2) as explained in Corollary 2.8
    beta12 = beta1*beta2;
    invsqbeta12 = (beta12*b0p*b0p + lc*(a0p-beta12)*(-a1p-beta1-beta2)*a0p*a0p)/(b0p*(b1p*a0p-a1p*b0p)*beta12+lc*((a0p-beta12)*a0p)**2);
    assert invsqbeta12**2 == 1/beta12, "sqrt is incorrect"

    # final transformation
    alpha = invsqbeta12;
    beta = -invsqbeta12*alpha2;
    gamma = 1;
    delta = -alpha1;
    epsilon =  (invsqbeta12*(alpha2-alpha1))**3;
    #auxiliary values to apply the transformation more efficiently
    al2 = alpha*alpha
    al2be = al2*beta
    al2de = al2*delta
    al3 = al2*alpha
    albe = alpha*beta
    alde = alpha*delta
    be2 = beta*beta
    bede = beta*delta
    de2 = delta*delta
    trafo = [alpha,beta,gamma,delta,epsilon,al2,al2be,al2de,al3,albe,alde,be2,bede,de2]


    Ap = (beta1+beta2)*invsqbeta12;
    Bp = (gamma1+gamma2)*invsqbeta12;
    Cp = gamma1*gamma2*invsqbeta12*invsqbeta12;
    Ep = 1/(lc*alpha);
    type1_invariants = [Ap,Bp,Cp,Ep]

    return type1_invariants, trafo


def apply_Mumford_transformation(P,trafo):
    """
    Transformation of point P in Mumford coordinates.
    On the level of curves, this cooresponds to
    x -> (alpha*x + beta)/(gamma*x + delta)
    y -> epsilon*y/(gamma*x + delta)**3
    Note: we assume gamma = 1
    And get some auxiliary information, e.g alde = alpha*delta
    """
    [al,be,ga,de,ep,al2,al2be,al2de,al3,albe,alde,be2,bede,de2] = trafo;    

    [a0,a1,a2,b0,b1] = P;

    if a2 == 0:
        print("not implemented");
    else:
        adeni = 1/(-a0 + de*a1 - de2)
        a1p = (al*(a0+a0) + (-alde - be)*a1 + (bede+bede))*adeni
        a0p = (-al2*a0 + albe*a1 - be2)*adeni
        b33 = -b0 + b1*de
        b22 =  al*(b0+b0+b0) - (be + alde+alde)*b1
        b11 = -al2*(b0+b0+b0) + (al2de + albe+albe)*b1
        b00 = al3*b0 - al2be*b1
        b1p = (a0p-a1p*a1p)*b33 + a1p*b22  - b11
        b0p = a0p*(-a1p*b33 + b22) - b00
    return [a0p,a1p,1,b0p,b1p]

def invariant_type2_type1(type2_invariants):
    #from type 2 invariants to type 1
    [A,B,C,E] = type2_invariants
    denom_aux = 1/((1-A)*(E-B+C))
    return [-(2+A+A)*(E-B+C)*denom_aux, (C+C-E-E)*(1-A)*denom_aux, (E+B+C)*(1-A)*denom_aux, (1-A)*(E-B+C)]

def invariant_type1_type2(type1_invariants):
    #from type 1 invariants to type 2
    [A,B,C,E] = type1_invariants
    term = -E*(A-2)
    return [(A+2)/(A-2), 2*term*(C-1),term*(B+C+1), term*(-B+C+1)]

def apply_trafo_type2_type1(P):
    #clearly, this special transformation can be evaluated more efficiently
    alpha = 1
    beta = -1
    gamma = 1
    delta = 1
    epsilon = 2
    trafo = [alpha,beta,gamma,delta,epsilon,1,-1,1,1,-1,1,1,-1,1]

    return apply_Mumford_transformation(P,trafo)

def apply_trafo_type1_type2(P):
    
    [a0,a1,a2,b0,b1] = P;

    if a2 == 0:
        print("not implemented");
    else:
        adeni = 1/(-a0 - a1 - 1)
        a1p = ((a0+a0)  - 2)*adeni
        a0p = (-a0 + a1 - 1)*adeni
        b33 = -b0 - b1
        b22 =  (b0+b0+b0) + b1
        b11 = - (b0+b0+b0) + b1
        b00 = b0 - b1
        b1p = (a0p-a1p*a1p)*b33 + a1p*b22  - b11
        b0p = a0p*(-a1p*b33 + b22) - b00
    return [a0p,a1p,1,b0p,b1p]


def type1_isogeny(type1_invariants):
    """
    compute the codomain invariants
    """
    [A,B,C,E] = type1_invariants
    if C == 1:
        print("codomain is product of elliptic curves.")
    denom_aux = 1/(1-C)
    Ap = C
    Bp = E + E
    Cp = (B-A*C)*E*denom_aux
    Ep = (A-B)*E*denom_aux

    return [Ap,Bp,Cp,Ep]

def auxiliary_products(type1_invariants):
    [A,B,C,E] = type1_invariants
    AB = A*B
    AC = A*C
    B2 = B*B
    BC = B*C
    C2 = C*C
    CE = C*E
    return [AB,AC,B2,BC,C2,CE]

def auxiliary_doubling(type2_invariants):
    [A,B,C,E] = type2_invariants
    t1 = (A+1)*E - C
    AC = A*C
    return [t1,AC]

def apply_type1_isogeny(P, type1_invariants,aux_iso):
    # Formula from Theorem 4.7

    [A,B,C,E] = type1_invariants
    [AB,AC,B2,BC,C2,CE] = aux_iso
    
    [a0,a1,a2,b0,b1] = P
    #computing monomials
    a02 = a0*a0
    a0a1 = a0*a1
    a12 = a1*a1
    a0b0 = a0*b0
    a0b1 = a0*b1
    a1b0 = a1*b0
    a12b0 = a12*b0
    a02a1b0 = a02*a1b0
    a0a12b0 = a0*a12b0
    a13b0 = a12*a1b0
    a02b1 = a02*b1
    a03b1 = a0*a02b1 
    a02a1b1 = a0a1*a0b1
    a0a12b1 = a0b1*a12
    a02b0 = a02*b0
    a0a1b0 = a0a1*b0
    a0a1b1 = a0a1*b1
    a021b1 = a02*b1
    b02 = b0*b0


    a00 = C2*(a02+a12) + BC*(a0a1+a1) + (B2 - C2 - C2)*a0 + C2
    a11 = (C2 - C)*(a0a1+a0a1+a1+a1) + (BC - B)*(a0+a0+a0+a0)
    a22 = -(C+C)*(a02+a12+1) - (BC + B)*(a0a1+a1) + (-B2 + C2 + C2 - C -C + 2)*(a0+a0)
    a33 = (1-C)*(a0a1+a0a1+a1+a1) + (B - BC)*(a0+a0+a0+a0)
    a44 = a02 + B*(a0a1+a1) + a12 + (B2 - 2)*a0 + 1
    b00 = C*(a02a1b0+a13b0-a03b1-a0a12b1- a0a1b0-a0a1b0+a02b1+a02b1+a1b0-a0b1) + B*(a0a12b0 - a02a1b1-a02b0+a0b0) + AC*(a02b0+a12b0-a0a1b1-a0b0) + AB*(a0a1b0-a02b1)      
    b11 = C*(a0a12b0+a0a12b0-a02a1b1-a02a1b1-a02b0+a12b0+b0) + (2)*(a02a1b1 - a0a12b0+a02b0-a0b0) + AC*(a0a1b0+a0a1b0-a02b1-a0b1+a1b0) + A*(-a0a1b0-a0a1b0+a02b1+a02b1) + B*(a0a1b0-a02b1+a0b1) + AB*a0b0
    b22 = - a02a1b0 - a13b0 + a03b1  + a0a12b1 - a1b0 - a0b1 + B*(-a0a12b0+a02a1b1+a02b0-a0b0) + A*(-a02b0-a12b0+a0a1b1-a0b0) + AB*(-a0a1b0+a02b1) + C*(a0a1b0+a0a1b0-a02b1-a02b1+a0b1+a0b1)  + AC*(a0b0+a0b0)
    b33 = - a02b0 - a12b0 + a0b0 + a0b0 - b0 + B*(-a0a1b0+a02b1-a0b1) + A*(-a02b1-a1b0+a0b1) - AB*a0b0
    bden =  (1 - a0) * (b02 - b1*(a1b0-a0b1));
    bden2 = bden*bden

    # Reduction Step
    lam = E*(B-A)*bden2 - (C-1)*b33*b33;
    #assert lam*a00 != 0, "edge case not implemented"
    
    a00lam = a00*lam
    den_aux = 1/(a00lam*bden)
    den = den_aux*bden
    a0p = (CE*bden2*(AC-B)+ (1-C)*b00*b00)*a44*den;
    a1p = (((2-C-C)*(CE*bden2 + b00*b11))*a44-a11*a0p*lam)*den;
    bdeni = den_aux*a00lam;
    b1p = (-b11+b22*a1p+b33*(a0p-a1p*a1p))*bdeni;
    b0p = (-b00+(b22-b33*a1p)*a0p)*bdeni;

    return [a0p,a1p,1,b0p,b1p]