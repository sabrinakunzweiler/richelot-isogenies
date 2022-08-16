SymplecticBasisg2 := function(J,m,cofactor)
  /*
  Symplectic Basis for the m-torsion of J (not deterministic)
  output: (P1,P2,Q1,Q2) with  e_m(P1,Q1) = e_m(P2,Q2) = zeta (primitive m-th root)
                              e_m(Pi,Qj) = 1 (for i neq j)
                              e_m(Pi,Pj) = e_m(Qi,Qj) = 1

  Assumptions: (m*cofactor)^2 = #J[FF_p],
   for example in the SIDH setting: m=2^ea, cofactor = 3^eb,
  */

  P:=[];
  Q:=[];
  repeat
      P[1] := cofactor * Random(J);
  until Order(P[1],2,m) eq m;

  repeat
  	P[3] := cofactor * Random(J);
    zeta := WeilPairing(P[1],P[3],m);
  until Order(zeta) eq m;
  Q[1] := P[3];

  repeat
    P[2] := cofactor * Random(J);
    alpha := [[Log(zeta,WeilPairing(P[i],P[j],m)): j in [1,2,3]]: i in [1,2,3]];
    P[2]:= P[2] - alpha[1][2]*P[3] + alpha[3][2]*P[1];
  until Order(P[2], 2, m) eq m;

  repeat
    P[4] := cofactor * Random(J);
    a34 := Log(zeta,WeilPairing(P[2],P[4],m));
  until GCD(m,a34) eq 1;

  t1,t2,t3 := XGCD(m,a34);
  P[4] := t3*P[4];

  for i in [1,2,3] do
    alpha[i][4] := Log(zeta,WeilPairing(P[i],P[4],m));
  end for;

  Q[2] := P[4] + alpha[3][4]*P[1] - alpha[1][4]*P[3];

  return [P[1],P[2],Q[1],Q[2]];
end function;

SymplecticBasisPowerTwo := function(J,power2,cofactor)
  /*
  Symplectic Basis for the m-torsion of J (not deterministic)
  output: (P1,P2,Q1,Q2) with  e_m(P1,Q1) = e_m(P2,Q2) = zeta (primitive m-th root)
                              e_m(Pi,Qj) = 1 (for i neq j)
                              e_m(Pi,Pj) = e_m(Qi,Qj) = 1

  Assumptions: (m*cofactor)^2 = #J[FF_p],
   for example in the SIDH setting: m=2^ea, cofactor = 3^eb,
  */
  m :=2^power2;
  P:=[];
  Q:=[];
  repeat
      P[1] := cofactor * Random(J);
  until 2^(power2-1)*P[1] ne Id(J);

  repeat
  	P[3] := cofactor * Random(J); //P[2]
    zeta := WeilPairing(P[1],P[3],m);
  until Order(zeta) eq m;
  Q[1] := P[3];

  repeat
    P[2] := cofactor * Random(J);
    alpha := [[Log(zeta,WeilPairing(P[i],P[j],m)): j in [1,2,3]]: i in [1,2,3]];
    P[2]:= P[2] - alpha[1][2]*P[3] + alpha[3][2]*P[1];
  until 2^(power2-1) *P[2] ne Id(J);

  repeat
    P[4] := cofactor * Random(J);
    a34 := Log(zeta,WeilPairing(P[2],P[4],m));
  until GCD(m,a34) eq 1;

  t1,t2,t3 := XGCD(m,a34);
  P[4] := t3*P[4];

  for i in [1,2,3] do
    alpha[i][4] := Log(zeta,WeilPairing(P[i],P[4],m));
  end for;

  Q[2] := P[4] + alpha[3][4]*P[1] - alpha[1][4]*P[3];
  return [P[1],P[2],Q[1],Q[2]];
end function;


SpecialBasisg2 := function(inverse_montgomery_invariants, power2, cofactor)
  /*
  input: Type-2 invariants defining hyperell. curve over Fp^2 with p = 2^power2 * cofactor -1
  output: (P1,P2,Q1,Q2) with  e_m(P1,Q1) = e_m(P2,Q2) = zeta (primitive m-th root)
                              e_m(Pi,Qj) = 1 (for i neq j)
                              e_m(Pi,Pj) = e_m(Qi,Qj) = 1

  Assumptions: (m*cofactor)^2 = #J[FF_p],
   for example in the SIDH setting: m=2^ea, cofactor = 3^eb,
  output: Symplectic Basis BB=(P1,P2,Q1,Q2) for J[2^n], s.t. 2^(n-1)*BB special 2-torsion basis.
  (ad-hoc implementation - not efficient)
  */
  m:= 2^power2;

  A := inverse_montgomery_invariants[1];
  B := inverse_montgomery_invariants[2];
  C := inverse_montgomery_invariants[3];
  d := inverse_montgomery_invariants[4];

  k := Parent(A);
  R<x> := PolynomialRing(k);
  Y := HyperellipticCurve((x^2-1)*(x^2-A)*(d*x^2-B*x+C));
  J := Jacobian(Y);

  BB := SymplecticBasisPowerTwo(J, power2,cofactor);
  BB2 := [2^(power2-1)*BB[1], 2^(power2-1)*BB[2], 2^(power2-1)*BB[3],2^(power2-1)*BB[4]];

  alpha := Sqrt(A);
  P1 := J ! [(x-1)*(x-alpha),0];
  Q1 := J ! [(x-1)*(x+1),0];
  if d eq 0 then
    beta := C/B;
    P2 := J ! [(x+alpha)*(x-beta),0];
    Q2 := J ! [(x-beta), 0];
  else
    delta := Sqrt(B^2-4*C*d);
    beta := (B + delta)/(2*d);
    gamma := (B-delta)/(2*d);
    P2 := J ! [(x+alpha)*(x-beta),0];
    Q2 := J ! [(x-beta)*(x-gamma), 0];
  end if;

  //base change from B2 to BB2 = (P1,P2,Q1,Q1)
  CC2 := [P1,P2,Q1,Q2];
  zeta := WeilPairing(CC2[1],CC2[3],2);
  check := Matrix(IntegerRing(), [[Log(zeta,WeilPairing(CC2[i],CC2[j],2)): j in [1,2,3,4]]: i in [1,2,3,4]]);
  A := Matrix(IntegerRing(), [[Log(zeta,WeilPairing(BB2[i],CC2[j],2)): j in [1,2,3,4]]: i in [3,4,1,2]]);
  //lift A from Z/2Z to Z/2^nZ

  B_new := [&+[A[i][j]*BB[i]: i in [1,2,3,4]]: j in [1,2,3,4]];
  //need to make B_new symplectic again
  zeta := WeilPairing(B_new[1], B_new[3],m);
  alpha := [[Log(zeta,WeilPairing(B_new[i],B_new[j],m)): j in [1,2,3]]: i in [1,2,3]];
  B_new[2] := B_new[2] + alpha[2][1]*B_new[3] + alpha[3][2]*B_new[1];
  a24 := Log(zeta,WeilPairing(B_new[2],B_new[4],m));
  t1,t2,t3 := XGCD(m,a24);
  B_new[4] := t3*B_new[4];
  for i in [1,2,3] do
    alpha[i][4] := Log(zeta,WeilPairing(B_new[i],B_new[4],m));
  end for;
  B_new[4] := B_new[4] + alpha[3][4]*B_new[1] - alpha[1][4]*B_new[3];

  return B_new;

end function;

RandomSymplecticGroup := function(B, power2)

  r1 := Random(2^power2-1);
  r2 := Random(2^power2-1);
  r3 := Random(2^power2-1);

  R1 := B[1] + r1*B[3] + r2*B[4];
  R2 := B[2] + r2*B[3] + r3*B[4];

  R := [R1,R2];
  R_coord := [];
  for P in R do
    a0 := Coefficient(P[1],0);
    a1 := Coefficient(P[1],1);
    a2 := Coefficient(P[1],2);
    b0 := Coefficient(P[2],0);
    b1 := Coefficient(P[2],1);

    R_coord := Append(R_coord, [a0,a1,a2,b0,b1]);
  end for;

  return R_coord;

end function;
