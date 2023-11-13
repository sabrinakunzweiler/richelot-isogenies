def find_special_trafo(invariants, kernel, point4):
    #we assume that kernel = <((x-alpha),0), x*(x-beta, 0)>,
    #point4 = [x0,x1,x2,x3] with 2*points4 = (x-alpha,0)
    [A,B,C,E] = invariants

    #avoid inversions
    large_denom = 1/(kernel[0][1]*kernel[1][0]*point4[0])
    alpha = kernel[0][2]*kernel[1][0]*point4[0]*large_denom
    beta = kernel[1][1]*kernel[0][1]*point4[0]*large_denom
    gamma = B - beta

    #symmteric polynomials and auxiliary symmetrics
    s2 = point4[1]*kernel[0][1]*kernel[1][0]*large_denom - alpha - alpha
    s4 = point4[2]*kernel[0][1]*kernel[1][0]*large_denom - (s2+alpha)*alpha
    t2 = 5*alpha - A - B
    t4 = 2*alpha*(3*alpha-B-B) - 3 + A*B + C
    prod = alpha*(alpha-beta)
    scalar = ((s4+(s2**2-t4)/2)*prod + (s4 - prod)*(s4 - prod))/(prod*(s2+s2+t2) + (s4 - prod)*(s2 + alpha + alpha - beta))
    #if not scalar.parent().is_field():
    #    scalar = lift_squareroot(prod, scalar, precision)
    assert scalar**2 == prod, "square-root in transformation is incorrect"
    scalar_inv = 1/scalar
    newA = (-2*alpha+beta)*scalar_inv #alphanew + 1/alphanew #
    newB = (A + gamma - 3*alpha)*scalar_inv #betanew + gammane
    newC = (A - 2*alpha)*(gamma - alpha)*scalar_inv*scalar_inv #betanew*gammanew
    newE = E*scalar

    #directly compute coefficients for the Kummer transformation
    a00 = prod
    a10 = -(alpha + alpha)*scalar  
    a11 = scalar
    a20 = alpha*alpha
    a21 = -alpha
    a30 = E*alpha*((alpha + alpha)*(-(alpha + alpha) + (A + B)) + (-A*B - C - 1))
    a31 = E*a20
    a32 = -E*(alpha +alpha)

    return [newA,newB,newC,newE], [a00,a10,a11,a20,a21,a30,a31,a32]

def apply_special_trafo(P, invariants, special_trafo):
    #transformation that sends the (2,2)-kernel into a kernel of the desired form.
    #It is defined by [a,b,c,d] = [1, -alpha, 0, scalar]
    [A,B,C,E] = invariants
    #[alpha,scalar] = special_trafo
    [a00,a10,a11,a20,a21,a30,a31,a32] = special_trafo
    
    [x0,x1,x2,x3] = P

    y0 = a00*x0
    y1 = a10*x0 + a11*x1
    y2 = a20*x0 + a21*x1 + x2
    y3 = a30*x0 + a31*x1 + a32*x2 + x3 
    #y0 = scalar**2*x0
    #y1 = scalar*(-2*alpha*x0 + x1)
    #y2 = alpha**2*x0 - alpha*x1 + x2
    #y3 = E*alpha*((2*alpha*(-2*alpha + (A + B)) + (-A*B - C - 1))*x0 + alpha*x1 - 2*x2) + x3
    
    return [y0,y1,y2,y3]


def invariant_type1_type2(type1_invariants):
    #from type 1 invariants to type 2
    [A,B,C,E] = type1_invariants
    denom = 1/(16*(A-2))
    term = (E*A - (E + E))*(A-2)*denom
    return [(A+2)*16*denom, 2*term*(1-C),term*(-B-C-1), term*(B-C-1)]

def compute_trafo_12(type1_invariants):
    [A,B,C,E] = type1_invariants
    a30 = -(A*C + B)*E
    a31 = A*B*E/2
    a32 = -E*(A+B)
    return [a30,a31,a32]

def apply_trafo_type1_type2(P, trafo_aux):
    [x0,x1,x2,x3] = P
    [x0D,x1D,x2D,x3D] = [x0+x0,x1+x1,x2+x2,x3+x3]
    [a30,a31,a32] = trafo_aux
    y0 = x0D - x1D + x2D
    y1 = x2D + x2D - x0D - x0D
    y2 = x0D + x1D + x2D
    y3 = a30*x0 + a31*x1 + a32*x2 + x3D
    #y3 = E*(-(A*C + B)*x0 +A*B/2*x1 - (A + B)*x2) + x3D
    return [y0,y1,y2,y3]

def Richelot_type2_invariants(type2_invariants):
    #return type1_invariants of the codomain
    [A,B,C,E] = type2_invariants
    Binv = 2/B
    newA = (E+C) * Binv
    newB = (A*E+C) * Binv
    newC = A
    newE = B + B
    return [newA,newB,newC,newE]

def Richelot_type1_invariants(type1_invariants):
     #return type2_invariants of the codomain
    [A,B,C,E] = type1_invariants
    deltainv = E/(4*(C-1))
    return [C, deltainv*(C+C-2) , deltainv*(A*C - B), deltainv*(B - A)]

def Richelot_type2_matrix(type2_invariants):
    #matrix defining the Richelot isogeny
    [A,B,C,E] = type2_invariants
    AB = A*B
    AC = A*C 
    AE = A*E
    B2 = B*B
    M0 = [AE - AC - C, 0  , 0 , 1 ,  C, -B,  0  ,  E, 0  ,0 ]
    M1 = [AB, -(AC + AC + AE + AE + C + C), AB + B, 0  , (AE + AE + C + C)*(C + E)/B, -(AE + AE + C + C + E + E),  1 ,  B, 0 , 0 ]
    M2 = [AC,  -AB,  0 , 0 , AE, 0  ,0  , -AE + C - E, 1  , 0 ]
    M3 = [A*(A*(4 *E*E -B2) - B2), 0 , 4 * A * (-B2  + (C + C) * E), 4*AE,  A*B2 , 0 , 0 , -B2*(A+1) + 4*C*C , C + C + C + C, 1 ]
    return [M0,M1,M2,M3]

def Richelot_type1_matrix(type1_invariants):
    #matrix defining the Richelot isogeny
    [A,B,C,E] = type1_invariants
    AC = A*C
    BC = B*C
    CC = C*C
    CE = C*E
    denom = 1/((CE + CE - E - E))
    delta = (CC - C - C + 1)*denom
    deltainv = E*E*denom
    t1 = AC - BC
    term = (AC + AC - BC  - B)*(AC + A - B - B) + (C + 1)*(CC - C - C +1) 
    L0 = [AC - BC, 0 , (-AC + B) + (-AC + B),  delta + delta, 0 ,  C - 1 , 0 , A - B, 0 , 0 ] # div 4
    L1 = [-CC + C, 0 , (-A*B + C + 1 )*(C-1), 0 , 0 , 0 , delta + delta, -C + 1, 0 , 0 ] #div 4
    L2 = [C*(-AC + B), CC - C, AC+AC-BC-BC,  0 ,  0 , 0 , 0 , -AC + B, delta +delta, 0 ] # div 4
    L3 = [deltainv*C*term, CE*(B - AC), -4 *deltainv*((AC-B)*(AC+AC-BC-BC) + C*(CC - C - C + 1)),  AC-BC , CE*(C-1)/2  , (AC-BC)*E, 0  ,  deltainv*term, B - AC, delta ]

    return [L0,L1,L2,L3]

def apply_Richelot_type2(P, M):
    
    [x0,x1,x2,x3] = P
    monomials = [x0*x0, x0*x1, x0*x2, x0*x3, x1*x1, x1*x2, x1*x3, x2*x2, x2*x3, x3*x3]
    y0 = M[0][0]*monomials[0] + M[0][3]*monomials[3] + M[0][4]*monomials[4] + M[0][5]*monomials[5] + M[0][7]*monomials[7]
    y1 = M[1][0]*monomials[0] + M[1][1]*monomials[1] + M[1][2]*monomials[2] + M[1][4]*monomials[4] + M[1][5]*monomials[5] + M[1][6]*monomials[6] + M[1][7]*monomials[7]
    y2 = M[2][0]*monomials[0] + M[2][1]*monomials[1] + M[2][4]*monomials[4] + M[2][7]*monomials[7] + M[2][8]*monomials[8]
    y3 = M[3][0]*monomials[0] + M[3][2]*monomials[2] + M[3][3]*monomials[3] + M[3][4]*monomials[4] + M[3][7]*monomials[7] + M[3][8]*monomials[8] + M[3][9]*monomials[9]
    
    return [y0,y1,y2,y3]

def apply_Richelot_type1(P, M):
    
    [x0,x1,x2,x3] = P
    monomials = [x0*x0, x0*x1, x0*x2, x0*x3, x1*x1, x1*x2, x1*x3, x2*x2, x2*x3, x3*x3]
    y0 = M[0][0]*monomials[0] + M[0][2]*monomials[2] + M[0][3]*monomials[3] + M[0][5]*monomials[5] + M[0][7]*monomials[7]
    y1 = M[1][0]*monomials[0] + M[1][2]*monomials[2] + M[1][6]*monomials[6] + M[1][7]*monomials[7]
    y2 = M[2][0]*monomials[0] + M[2][1]*monomials[1] + M[2][2]*monomials[2] + M[2][7]*monomials[7] + M[2][8]*monomials[8]
    y3 = M[3][0]*monomials[0] + M[3][1]*monomials[1] + M[3][2]*monomials[2] + M[3][3]*monomials[3] + M[3][4]*monomials[4] + M[3][5]*monomials[5] + M[3][7]*monomials[7] + M[3][8]*monomials[8] + M[3][9]*monomials[9]

    return [y0+y0,y1+y1,y2+y2,y3]

