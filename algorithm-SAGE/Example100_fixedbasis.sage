"""
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

In this example we compute a  (2^51,2^51)-isogeny starting at a predefined abelian surface.
This abelian surface is the Jacobian of a hyperelliptic curve defined by a Type-2 equation.
The (2^51,2^51)-group defining this isogeny is randomly selected.
"""

load("setup.sage");
load("Richelot_formulae.sage");

# typical G2SIDH parameters  with a 100-bit prime
p = 2^51*3^32-1;
Fp = FiniteField(p);
R.<x> = PolynomialRing(Fp);
Fp2.<z2> = Fp.extension(x^2+863349905383196036096441154452*x + 7);

type2_invariants = [
    2506647832581528277950979545255*z2 + 1361205151906637044440651372441,
    3342910751114468500993201630164*z2 + 3497607818217780653905550638544,
    1077077857247676286359122984243*z2 + 104569387409790195049388192267,
    302754864461984374458169335566*z2 + 219814990121629982473457168231
];

[A,B,C,E] = type2_invariants;

R.<x> = PolynomialRing(Fp2);
Y = HyperellipticCurve((x^2-1)*(x^2-A)*(E*x^2-B*x+C));
J = Jacobian(Y);

"""
A special symplectic basis was precomputed to speed-up the computation.
Alternatively it can be generated using the functions in "setup.m".
"""

BB = [ J(x^2 + (3220565603736012661433518932215*z2 +
    2153643323630220719758710398803)*x + 3037024312445676383637796554558*z2 +
    56474981936235249709409175773, (2755940572414405269664673908776*z2 +
    3136196238971113020834601483132)*x + 1090037494367476681437808853457*z2 +
    4117364898884110809454223452783),
    J(x^2 +
    (1264312533772572565211084233155*z2 + 491418407023645668000226288462)*x +
    2706673070197391863094470415341*z2 + 2359221575715870192970117021508,
    (3962267309606171003612387423816*z2 + 4171884226675813798823669787498)*x +
    2140839758022537943976774616930*z2 + 1448941331592607266475027054061),
    J(x^2 + (2158615300062249032631999429778*z2 + 2829690642397134608443113540492)*x
    + 1403134488571232131810296779726*z2 + 470479498660567317727613683176,
    (1185824149593457212262241895533*z2 + 4172499238468688860142878587861)*x +
    254860736533019469370581247850*z2 + 2314793371746776420544286643052),
    J(x^2 + (628395930236594467800939384709*z2 + 2952804811330269577546599691065)*x +
    614560107981011505766062756506*z2 + 2623100418911082082974286667437,
    (3172708646303957902517741449372*z2 + 3546134881676385321176106285451)*x +
    532101508677234858400718097420*z2 + 2183512758630963371074851515760)];

"""
In order to simulate the first round of a G2SIDH key exchange, where Alice needs
to compute the images of the 3^m-torsion, we generate 4 random points on the Jacobian.
"""

PP = [];
def random_point_jac(J, order=None):
    f = J.curve().hyperelliptic_polynomials()[0]
    counter = 0 #to make the choice of the sign of the square roots random
    while True:
        rt1 = Fp.random_element()
        if f(rt1).is_square():
            break
        counter += 1
    y1 = (-1)^counter*sqrt(f(rt1))
    while True:
        rt2 = Fp.random_element()
        if f(rt2).is_square():
            break
        counter += 1
    y2 = (-1)^counter*sqrt(f(rt2))
    a = (x-rt1)*(x-rt2)
    b = (y1-y2)/(rt1-rt2)*(x-rt1) + y1
    return J([a,b])

for i in range(4):
  P = random_point_jac(J);
  [a0,a1,a2] = P[0].coefficients()
  [b0,b1] = P[1].coefficients()
  PP.append([a0,a1,a2,b0,b1]);

ker = RandomSymplecticGroup(BB, 51);
t1 = cputime();
Ynew = RichelotChain(type2_invariants, ker, 51, P_list=PP, partition=[10,8,7,6,5,5], last_step=true);
print("CPU time for Richelot chain (Round 1):", cputime(t1))
t2 = cputime();
Ynew = RichelotChain(type2_invariants, ker, 51, partition=[10,8,7,6,5,5], last_step=true);
print("CPU time for Richelot chain (Round 2):", cputime(t2))
