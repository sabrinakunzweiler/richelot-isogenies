"""
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

In this example we compute a  (2^51,2^51)-isogeny starting at a predefined abelian surface.
This abelian surface is the Jacobian of a hyperelliptic curve defined by a Type-2 equation.
The (2^51,2^51)-group defining this isogeny is randomly selected.
"""

load("setup.sage");
load("Richelot_formulae.sage");
load("richelot_aux_PopeOudompheng.sage")

# typical G2SIDH parameters  with a 100-bit prime
p = 2^51*3^32-1;
#Fp = FiniteField(p);
Fp2.<z2> = GF(p^2, modulus=x^2+863349905383196036096441154452*x + 7)
#Fp2 = FiniteField(p^2);
#R.<x> = PolynomialRing(Fp);
#Fp2.<z2> = Fp.extension(x^2+863349905383196036096441154452*x + 7);
#z2 = Fp2.gen()
assert -z2^2 == 863349905383196036096441154452*z2 + 7
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
to compute the images of the 3^m-torsion, we generate 4 points on the Jacobian.
"""

PP = [];
def random_point_jac(J, order=None):
    #note: this only outputs point with K-rational support!
    f = J.curve().hyperelliptic_polynomials()[0]
    counter = 0 #to make the choice of the sign of the square roots random
    while True:
        rt1 = Fp2.random_element()
        if f(rt1).is_square():
            break
        counter += 1
    y1 = (-1)^counter*sqrt(f(rt1))
    while True:
        rt2 = Fp2.random_element()
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

print("Our algorithm:")
t1 = cputime();
inv_new, PP_new = RichelotChain(type2_invariants, ker, 51, P_list=PP, partition=[10,8,7,6,5,5], last_step=true);
print("CPU time for Richelot chain (Round 1):", cputime(t1))
t2 = cputime();
inv_new, _ = RichelotChain(type2_invariants, ker, 51, partition=[10,8,7,6,5,5], last_step=true);
print("CPU time for Richelot chain (Round 2):", cputime(t2))
[Ap,Bp,Cp,Ep] = inv_new;
Ynew = HyperellipticCurve((x^2-1)*(x^2-Ap)*(Ep*x^2-Bp*x+Cp));

print("The Pope-Oudompheng algorithm")
h = (x^2-1)*(x^2-A)*(E*x^2-B*x+C)
a = 51
kernel = [ker[0][0]+ker[0][1]*x+x^2,ker[0][3]+ker[0][4]*x,ker[1][0]+ker[1][1]*x+x^2,ker[1][3]+ker[1][4]*x];
points = [];
for P in PP:
    [a0,a1,a2,b0,b1] = P
    D = [x^2+a1*x+a0, b1*x+b0]
    points.append(D)

t1 = cputime();
hnew = IsogenyChainPO(h,kernel,points,a);
print("CPU time for Richelot chain (Round 1): ", cputime(t1))
t2 = cputime();
hnew = IsogenyChainPO(h,kernel,[],a);
print("CPU time for Richelot chain (Round 2): ", cputime(t2))
Ynew2 = HyperellipticCurve(hnew[0]);

print("  ");
print("171-bit prime");

p = 2^87*3^49*5^3-1;
Fp = FiniteField(p);
R.<x> = PolynomialRing(Fp);
Fp2.<z2> = Fp.extension(x^2+1)

type2_invariants =  [
    3446798252521358094517040152788188703411472634458348*z2 +
        3934251785801114798773063613093606362278332002469273,
    14442707041689001196121281746286164295977084155961*z2 +
        2109612034623761207119749551680295829371348991735167,
    812627276298732058211077115034347230811721036500113*z2 +
        772992546050249836417619075409507366069466860817407,
    4149184472552449328486684911846181278057937206823244*z2 +
        76217450280973533167781361110947786357072578100390
];

[A,B,C,E] = type2_invariants;

R.<x> = PolynomialRing(Fp2);
Y = HyperellipticCurve((x^2-1)*(x^2-A)*(E*x^2-B*x+C));
J = Jacobian(Y);


#A special symplectic basis was precomputed to speed-up the computation.
#Alternatively it can be generated using the functions in "setup.m".


BB = [ J(x^2 + (4228467102363626818406821291187317314387242721137631*z2 +
    3804540158078856260408798951914626839615188809113971)*x +
    4275797638838160436292168115677773824418883521347906*z2 +
    1137321469212884956492608367539267130850696607667781,
    (4610281969133590881952445078108252720655017745022532*z2 +
    3973686199211383520298404807386383773962007438550548)*x +
    648258542931357692710711566778124079589290457474914*z2 +
    215955646752530011869381596473856112953124938376677),
    J(x^2 +
    (4338911201117717228685979630682819315648743416178581*z2 +
    164719356863546810047096688631436661008130416921123)*x +
    682832746418814828949968317152011216794405052776050*z2 +
    3227347167670773728077792441891694104440549509869116,
    (3338502728879829955316606942134401994124769140969765*z2 +
    4592275636898486852119197503541440469250199985025556)*x +
    2757504835057898567366283741674487382402327863649680*z2 +
    3009521334788140091333540237394896005087239208970452),
    J(x^2 +
    (2555738952316265353748094437339089294680814610537684*z2 +
    2176640765922335894797016044059634690560411512928873)*x +
    2458179105692663378496708693061286763272601074742808*z2 +
    1624818658770637954475518067668700878207038193447942,
    (4088437804309543062900709175013586494723123742583908*z2 +
    2807979495610579340461439943471276151550973301745513)*x +
    955764076954805838607110282177328845041233230865666*z2 +
    564011195494660755723258516073314725149407100733060),
    J(x^2 +
    (3764355748013851996734964908004966244986492767609594*z2 +
    4552431470848504784464885569352125742862744106050259)*x +
    3689157838657674380945494042451165029342922168026393*z2 +
    2127491357759698556236561733691614246364243090239406,
    (916752860733121206178059110704899351395203541469843*z2 +
    1183439552501690970580496543748515900250690676492295)*x +
    128074129624897914327260173767029947581651278749621*z2 +
    1819697913691524214918490466769661929366215392418827) ];

PP = []
for i in range(4):
  P = random_point_jac(J);
  [a0,a1,a2] = P[0].coefficients()
  [b0,b1] = P[1].coefficients()
  PP.append([a0,a1,a2,b0,b1]);

ker = RandomSymplecticGroup(BB,87); #slow!

print("Our algorithm")
t1 = cputime();
Ynew = RichelotChain(type2_invariants, ker,87,P_list=PP,partition=[12,11,10,9,8,8,7,6,6], last_step=true);
print("CPU time for Richelot chain (Round 1): ",  cputime(t1))
t2 = cputime();
Ynew = RichelotChain(type2_invariants, ker,87,partition=[12,11,10,9,8,8,7,6,6],last_step=true);
print("CPU time for Richelot chain (Round 2): ",  cputime(t2))


print("The Pope-Oudompheng algorithm")
h = (x^2-1)*(x^2-A)*(E*x^2-B*x+C)
a = 87
kernel = [ker[0][0]+ker[0][1]*x+x^2,ker[0][3]+ker[0][4]*x,ker[1][0]+ker[1][1]*x+x^2,ker[1][3]+ker[1][4]*x];
points = [];
for P in PP:
    [a0,a1,a2,b0,b1] = P
    D = [x^2+a1*x+a0, b1*x+b0]
    points.append(D)

t1 = cputime();
hnew = IsogenyChainPO(h,kernel,points,a);
print("CPU time for Richelot chain (Round 1): ", cputime(t1))
t2 = cputime();
hnew = IsogenyChainPO(h,kernel,[],a);
print("CPU time for Richelot chain (Round 2): ", cputime(t2))
Ynew2 = HyperellipticCurve(hnew[0]);
