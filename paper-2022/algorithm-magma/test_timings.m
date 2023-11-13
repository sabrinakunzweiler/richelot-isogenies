/*
AUTHOR: Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>

In this example we compute a  (2^51,2^51)-isogeny starting at a predefined abelian surface.
This abelian surface is the Jacobian of a hyperelliptic curve defined by a Type-2 equation.
The (2^51,2^51)-group defining this isogeny is randomly selected.
 */

load "setup.m";
load "Richelot_formulae.m";
load "richelot_aux_CastryckDecru.m";


/* typical G2SIDH parameters  with a 100-bit prime*/
print "100-bit prime";
p := 2^51*3^32-1;
Fp := FiniteField(p);
Fp2<i> := FiniteField(p^2);


type2_invariants := [
    864910784845209984294141672783*i + 1841246881133239550833588853065,
    1531260130020582574540405700414*i + 839920513885402936304320140616,
    1577799312233233169016154406767*i + 3876628765285855606468095639872,
    3347647507247385936903515809427*i + 283879061111347919579893365565
];

A,B,C,E := Explode(type2_invariants);

R<x> := PolynomialRing(Fp2);
Y := HyperellipticCurve((x^2-1)*(x^2-A)*(E*x^2-B*x+C));
J := Jacobian(Y);

/*
A special symplectic basis was precomputed to speed-up the computation.
Alternatively it can be generated using the functions in "setup.m".
*/
BB := [ J! [x^2 + (2048177950127403145338518498193*i + 3260855624648414854597815770936)*x
    + 1065943272518777363189308359849*i + 4170059877096463464170504064042,
    (1787473999945136415372677611941*i + 3225174621324581946500416415121)*x +
    1053736896742465439529080350931*i + 1161769639239578895072699582562],
    J! [x^2 + (3262648533685753908205101802985*i + 2412333043753863750404843555436)*x +
    3299615980789628215290144544695*i + 2357008346492464855870220668355,
    (3050658048693290942615297001459*i + 565809300240317436798550170525)*x +
    1115620414697258498624766818695*i + 3072981303424689146560429934233],
    J! [x^2 + (1433432644451892360855296005961*i + 4111175860429432028586250312861)*x +
    3650309464994401871670899704281*i + 1396238902566833541466344899335,
    (3347745822015355836307248310637*i + 2459320396780179782488149760336)*x +
    2066013158708777508116234304814*i + 2611357719910616559217955288830],
    J! [x^2 + (4025901915813799414477572929381*i + 1372764163987217502262279361663)*x +
    4150191959457132531287089077796*i + 1161371131629838050908038089038,
    (869768689127460178856713756976*i + 2994839532602888699614596271321)*x +
    2229405191288023738261851781672*i + 3697829497504888106334393713911] ];

/*
In order to simulate the first round of a G2SIDH key exchange, where Alice needs
to compute the images of the 3^m-torsion, we generate 4 random points on the Jacobian.
*/

PP := [];
for i := 1 to 4 do
  P := Random(J);
  a0 := Coefficient(P[1],0);
  a1 := Coefficient(P[1],1);
  a2 := Coefficient(P[1],2);
  b0 := Coefficient(P[2],0);
  b1 := Coefficient(P[2],1);
  PP := Append(PP, [a0,a1,a2,b0,b1]);
end for;

ker := RandomSymplecticGroup(BB, 51);

t1 := Cputime();
Ynew := RichelotChain(type2_invariants, ker,51:P_list:=PP);
printf "CPU time for Richelot chain (Round 1): %o\n",  Cputime(t1);
t2 := Cputime();
Ynew := RichelotChain(type2_invariants, ker,51);
printf "CPU time for Richelot chain (Round 2): %o\n",  Cputime(t2);
print "Timing for the Castryck-Decru algorithm:";


print "Timing for the Castryck-Decru algorithm:";
//translation to the notation in the C-D algorithm
A,B,C,E := Explode(type2_invariants);
h := (x^2-1)*(x^2-A)*(E*x^2-B*x+C);
a := 51;
kernel := [ker[1][1]+ker[1][2]*x+x^2,ker[1][4]+ker[1][5]*x,ker[2][1]+ker[2][2]*x+x^2,ker[2][4]+ker[2][5]*x];
points := [];
for P in PP do
    a0,a1,a2,b0,b1 := Explode(P);
    D := [x^2+a1*x+a0, b1*x+b0];
    points := Append(points, D);
end for;

t1 := Cputime();
hnew := IsogenyChainCD(kernel,points,h,a);
printf "CPU time for Richelot chain (Round 1): %o\n", Cputime(t1);
t2 := Cputime();
hnew := IsogenyChainCD(kernel,[],h,a);
printf "CPU time for Richelot chain (Round 2): %o\n", Cputime(t2);

print "  ";
print "171-bit prime";
/* typical G2SIDH parameters  with a 171-bit prime (expected 128-bit security)*/
p := 2^87*3^49*5^3-1;
Fp2<z2> := FiniteField(p^2);

type2_invariants :=  [
    3446798252521358094517040152788188703411472634458348*z2 +
        3934251785801114798773063613093606362278332002469273,
    14442707041689001196121281746286164295977084155961*z2 +
        2109612034623761207119749551680295829371348991735167,
    812627276298732058211077115034347230811721036500113*z2 +
        772992546050249836417619075409507366069466860817407,
    4149184472552449328486684911846181278057937206823244*z2 +
        76217450280973533167781361110947786357072578100390
];

A,B,C,E := Explode(type2_invariants);

R<x> := PolynomialRing(Fp2);
Y := HyperellipticCurve((x^2-1)*(x^2-A)*(E*x^2-B*x+C));
J := Jacobian(Y);

/*
A special symplectic basis was precomputed to speed-up the computation.
Alternatively it can be generated using the functions in "setup.m".
*/

BB := [ J![x^2 + (4228467102363626818406821291187317314387242721137631*z2 +
    3804540158078856260408798951914626839615188809113971)*x +
    4275797638838160436292168115677773824418883521347906*z2 +
    1137321469212884956492608367539267130850696607667781,
    (4610281969133590881952445078108252720655017745022532*z2 +
    3973686199211383520298404807386383773962007438550548)*x +
    648258542931357692710711566778124079589290457474914*z2 +
    215955646752530011869381596473856112953124938376677],
    J! [x^2 +
    (4338911201117717228685979630682819315648743416178581*z2 +
    164719356863546810047096688631436661008130416921123)*x +
    682832746418814828949968317152011216794405052776050*z2 +
    3227347167670773728077792441891694104440549509869116,
    (3338502728879829955316606942134401994124769140969765*z2 +
    4592275636898486852119197503541440469250199985025556)*x +
    2757504835057898567366283741674487382402327863649680*z2 +
    3009521334788140091333540237394896005087239208970452],
    J![x^2 +
    (2555738952316265353748094437339089294680814610537684*z2 +
    2176640765922335894797016044059634690560411512928873)*x +
    2458179105692663378496708693061286763272601074742808*z2 +
    1624818658770637954475518067668700878207038193447942,
    (4088437804309543062900709175013586494723123742583908*z2 +
    2807979495610579340461439943471276151550973301745513)*x +
    955764076954805838607110282177328845041233230865666*z2 +
    564011195494660755723258516073314725149407100733060],
    J![x^2 +
    (3764355748013851996734964908004966244986492767609594*z2 +
    4552431470848504784464885569352125742862744106050259)*x +
    3689157838657674380945494042451165029342922168026393*z2 +
    2127491357759698556236561733691614246364243090239406,
    (916752860733121206178059110704899351395203541469843*z2 +
    1183439552501690970580496543748515900250690676492295)*x +
    128074129624897914327260173767029947581651278749621*z2 +
    1819697913691524214918490466769661929366215392418827] ];

/*
In order to simulate the first round of a G2SIDH key exchange, where Alice needs
to compute the images of the 3^m-torsion, we generate 4 random points on the Jacobian.
*/

PP := [];
for i := 1 to 4 do
      P := Random(J);
      a0 := Coefficient(P[1],0);
      a1 := Coefficient(P[1],1);
      a2 := Coefficient(P[1],2);
      b0 := Coefficient(P[2],0);
      b1 := Coefficient(P[2],1);
      PP := Append(PP, [a0,a1,a2,b0,b1]);
end for;

ker := RandomSymplecticGroup(BB,87);


t1 := Cputime();
Ynew := RichelotChain(type2_invariants, ker,87:P_list:=PP);
printf "CPU time for Richelot chain (Round 1): %o\n",  Cputime(t1);
t2 := Cputime();
Ynew := RichelotChain(type2_invariants, ker,87);
printf "CPU time for Richelot chain (Round 2): %o\n",  Cputime(t2);

print "Timing for the Castryck-Decru algorithm:";
//translation to the notation in the C-D algorithm
A,B,C,E := Explode(type2_invariants);
h := (x^2-1)*(x^2-A)*(E*x^2-B*x+C);
a := 87;
kernel := [ker[1][1]+ker[1][2]*x+x^2,ker[1][4]+ker[1][5]*x,ker[2][1]+ker[2][2]*x+x^2,ker[2][4]+ker[2][5]*x];
points := [];
for P in PP do
    a0,a1,a2,b0,b1 := Explode(P);
    D := [x^2+a1*x+a0, b1*x+b0];
    points := Append(points, D);
end for;

t1 := Cputime();
hnew := IsogenyChainCD(kernel,points,h,a);
printf "CPU time for Richelot chain (Round 1): %o\n", Cputime(t1);
t2 := Cputime();
hnew := IsogenyChainCD(kernel,[],h,a);
printf "CPU time for Richelot chain (Round 2): %o\n", Cputime(t2);
