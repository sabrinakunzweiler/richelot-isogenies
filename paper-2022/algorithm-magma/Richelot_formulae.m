/*
* Code accompanying the paper "Efficient Computation of (2^n,2^n)-Isogenies"
* AUTHOR: Sabrina Kunzweiler
*
* (C) 2022 Sabrina Kunzweiler <sabrina.kunzweiler@ruhr-uni-bochum.de>
* Distributed under the terms of the GNU General Public License
* https://www.gnu.org/licenses/*
*/


DoubleType2 := function(P_list, type2_invariants)
  /*
  * INPUT:  -   type2_invariants = [A,B,C,E] that define the genus-2 curve
  *              y^2 = (x^2-1)(x^2-A)(Ex^2-Bx+C)
  *         -   P_list: list of points on the Jacobian of the genus-2 curve
  * OUTPUT: a list of points containing the doubles of the points in P_list
  */

  A,B,C,E := Explode(type2_invariants);

  /* Precomputations */
  t1 := (A+1)*E - C;
  AC := A*C;

  P_list_new := [];
  for P in P_list do
    a0,a1,a2,b0,b1 := Explode(P);
    if a2 eq 0 then
      print "not implemented";
    end if;
    /* Precomputations */
    a1b1 := a1*b1;
    a02:= a0^2;
    a12 := a1^2;
    b02 := b0^2;
    a0b1 := a0*b1;
    a1b0 := a1*b0;

    /* Composition step: 2*P = (a^2, bb), where bb = b + a * (e1x+e0) */
    Edeni := 1/(2*a0*(b1*(a1b0 - a0b1) - b02));
    E0 := a02*(E*(a12*(a1b1 + 3*b0) -2*a0*(2*a1b1 + b0))+ B*((a12-A-1)*b1  + 2*(a1b0 -a0b1)) -t1*(a1b1+b0)) + (AC-b02)*(a1b1-b0);
    E1 := a02*( 2*E*(-3*a1b0  + a0b1) - (3*E*a1 +2*B)*a1b1)  + a0*(4*E*a12*a1b0 +t1*(a0b1-2*a1b0) - (2*a0+A+1-3*a12)*B*b0) + b1*(AC- b02);
    /* formulae above only work for a0 != 0, the formulae below are more general but slightly slower */
    //Edeni :=  1/(2*(b1*(a0b1-a1b0) + b02));
    //E1 := E*(a12*(-a12*b1 - 4*a1b0 + 6*a0b1) + 3*a0*(2*a1b0 - a0b1))+ B*(a1*(-a12*b1 -3*a1b0 + 4*a0b1) +a0b02)  + t1*(a12*b1 +2*(a1b0-a0b1)) + b1*(b12 -t2)  + (AB + B)*(a1b1+b0);
    //E0 := E*(a12*(-3*a0b02 + (a1b0 - 2*a0b1)*a1)+ 3*a02*(2*a1b1  + b0)) + B*(a1b0*(a12 -4*a0 -1) +a0b1*(-2*a12+3*a0+2)) + t1*(a1*(2*a0b1 -a1b0) + a0b02) + b0*(b12 +t2)  + AB*(2*a0b1-a1b0 +b1);
    e0 := E0*Edeni;
    e1 := E1*Edeni;
    b00 := a0*e0+b0;
    b11 := a1*e0+a0*e1+b1;
    b22 := a1*e1+e0;
    b33 := e1;

    /* Reduction Step */
    b332 := b33^2;
    den := E - b332;
    if den eq 0 then
      a0p := (-AC + b00^2)/(a02*(B + 2*b22*b33));
      a1p := 1;
      a2p := 0;
      b0p := (b33*a0p - b22)*a0p^2 + b11*a0p - b00;
      b1p := 0;
    else
      a0p := (AC - b00^2)/(a02*den);
      a1p := (2*a1*(b332-E) - B  - 2*b22*b33)/den;
      a2p := 1;
      b1p :=  (-b33*(a1p^2 - a0p) + b22*a1p - b11);
      b0p := (-b33*a0p*a1p + b22*a0p - b00);
    end if;

    P_list_new := Append(P_list_new,[a0p,a1p,a2p,b0p,b1p]);
  end for;
  return P_list_new;

end function;


MumfordTransformation := function(trafo, P_points)
  al,be,ga,de,ep := Explode(trafo);
  Qpoints := [];

  alde := al*de;
  alga := al*ga;
  bega := be*ga;
  bede := be*de;

  for P in P_points do
    a0,a1,a2,b0,b1 := Explode(P);

    if a2 eq 0 then
      print "not implemented";
    else
      adeni := 1/((-a0*ga + a1*de)*ga - de^2);
      a1p := (2*(a0*alga+bede) - a1*(bega + alde))*adeni;
      a0p := ((-a0*al + a1*be)*al - be^2)*adeni;
      b33 := (-b0*ga + b1*de)*ga^2;
      b22 :=  (3*b0*alga - b1*(bega + 2*alde))*ga;
      b11 := (-3*b0*alga + b1*(2*bega + alde))*al;
      b00 := (b0*al - b1*be)*al^2;
      b1p := (a0p-a1p^2)*b33 + a1p*b22  - b11;
      b0p := a0p*(-a1p*b33 + b22) - b00;
    end if;

    Qpoints := Append(Qpoints, [a0p,a1p,1,b0p,b1p]);
  end for;

  return Qpoints;
end function;

MumfordTransformationga1 := function(trafo, P_lists)
  //in our applications ga =1
  al,be,ga,de,ep := Explode(trafo);
  alde := al*de;
  bede := be*de;

  Q_lists := [];
  for P_points in P_lists do
    Q_points := [];
    for P in P_points do
      a0,a1,a2,b0,b1 := Explode(P);

      if a2 eq 0 then
        print "not implemented";
      else
        adeni := 1/(-a0 + (a1-de)*de);
        a1p := (2*(a0*al+bede) - a1*(be + alde))*adeni;
        a0p := ((-a0*al + a1*be)*al - be^2)*adeni;
        b33 := -b0 + b1*de;
        b22 :=  3*b0*al - b1*(be + 2*alde);
        b11 := (-3*b0*al + b1*(2*be + alde))*al;
        b00 := (b0*al - b1*be)*al^2;
        b1p := (a0p-a1p^2)*b33 + a1p*b22  - b11;
        b0p := a0p*(-a1p*b33 + b22) - b00;
      end if;
      Q_points := Append(Q_points, [a0p,a1p,1,b0p,b1p]);
    end for;
    Q_lists := Append(Q_lists, Q_points);
  end for;

  return Q_lists;
end function;



RichelotType1 := function(type1_invariants, P_lists)
  /*
  * INPUT:  - type1_invariants = [A,B,C,E]: defining a genus-2 curve
  *           y^2 = E x (x^2-A*x+1)(x^2-B*x+C)
  *        -  P_lists: lists of lists of points on the jacobian of the genus-2 curve
  * OUTPUT: - type2_invariants = [A',B',C',E'] defining
  *            y^2 = (x^2-1)*(x^2-A')*(E'x^2-B'*x+C')
  *           such that its Jacobian is the codomain of the isogeny phi with kernel
  *           G = <(x,0),(x^2-A*x+1)>
  *         - Q_list: lists of points s.t. phi(P)=Q for any point P in P_lists

  * Assumptions: C != 1  + assumptions from Theorem 4.7.
  */

  A, B, C, E := Explode(type1_invariants);
  if C eq 1 then
    print "not implemented: codomain is product of elliptic curves.";
  end if;

  Q_lists := [];
  for P_list in P_lists do
    Q_list := [];
    for P in P_list do
      a0,a1,a2,b0,b1 := Explode(P);
      if a2 eq 0 then
        print "not implemented";
      end if;

      /* Formula from Theorem 4.7*/
      mu := a1*b0-a0*b1;
      t1 := (a0-1)^2+a1^2;
      t2 := (a0+1)*a1;
      a00 :=  C*(C*t1 + B*t2) +  B^2*a0;
      a11 :=  (C - 1)*(2*C*t2 + 4*B*a0);
      a22 :=  -2*C*t1 - B*(C +1)*t2 + 2*(-B^2 + 2*(C-1)^2)*a0;
      a33 :=  2*(-C + 1)*(t2 +2*B*a0);
      a44 :=  t1 + B*(t2 +  B*a0);
      b00 :=  mu*(C*(t1+A*a1) + B*a0*(a1+A)) + a0*b0*(a0-1)*(A*C-B);
      b11 :=  a0*(2*(C - 1)*(a1*mu-a0*b0)  + ((C- 2)*A+B)*mu +  (A*B - 2)*b0 + B*b1) +C*(A*((t2-a1)*b0 + mu) +b0*(t1+2*a0));
      b22 :=  -1*(t1  + B*(t2-a1)  + A*(B*a0 + a1))*mu+ a0*((B-A)*(a0-1)*b0  + 2*(C-1)*(b0*A+b1+mu));
      b33 :=  -( (A*a0+t2)*B+t1)*b0  + (B-A)*(mu+a0^2*b1);
      bden :=   (1 - a0) * (-b1*mu + b0^2);

      /* Reduction Step */
      lam := E*(B-A)*bden^2 - (C-1)*b33^2;
      if lam*a00 eq 0 then
        print "not implemented";
        Q_list := Append(Q_list, [0,0,0,0,0]);
      else
        den := 1/(lam*a00);
        a0p := (E*bden^2*C*(A*C-B)+ (1-C)*b00^2)*a44*den;
        a1p := (((1-C)*(E*bden^2*2*C + 2*b00*b11))*a44-a11*a0p*lam)*den;
        bdeni := 1/bden;
        b1p := (-b11+b22*a1p+b33*(a0p-a1p^2))*bdeni;
        b0p := (-b00+(b22-b33*a1p)*a0p)*bdeni;

        Q_list := Append(Q_list, [a0p,a1p,1,b0p,b1p]);
      end if;
    end for;
    Q_lists := Append(Q_lists, Q_list);
  end for;
  Ap := C;
  Bp := 2*E;
  Cp := (B-A*C)*E/(1-C);
  Ep := (A-B)*E/(1-C);

  return [Ap,Bp,Cp,Ep], Q_lists;

end function;


FindTransformation := function(type2_invariants, g1, g2:  P:=[])
  /*
  * INPUT:  - type2_invariants =[A,B,C,E] defining the hyperelliptic curve
  *           y^2 = (x^2-1)(x^2-A)(Ex^2-Bx+C)
  *         - g1 =[g10,g11,g12], g2 = [g20,g21,g22]
  *           coefficients of quadratic polynomials g1,g2 satisfying:
  *               - g1*g2 divide (x^2-1)(x^2-A)(Ex^2-Bx+C),
  *               - g1,g2 are not equal to any of the quadratic factors,
  *               - g1(1)=0 or g1(-1)=0
  *         - P = [a0,a1,a2,b0,b1] a point of order 4 with 2*P = (g1,0).
              this is optional (only used to speedup the computations)
  * OUTPUT: - [alpha,beta,gamma,delta,epsilon] defining the coordinate transformation
  *            from Proposition 4.1 to a Type-1 equation
  *         - [Ap,Bp,Cp,Ep] coeffcients of the Type-1 equation
  */

  A,B,C,E := Explode(type2_invariants);
  g10,g11,g12 := Explode(g1);
  g20,g21,g22 := Explode(g2);
  /* proceed as in Proof of Proposition 4.5 */
  if g10 eq -g11-1 then //test if 1 is a root of g1
    alpha1 := 1;
    alpha2 := g10;
  else //-1 is a root
    alpha1 := -1;
    alpha2 := -g10;
  end if;

  /*
  * either -alpha1 or -alpha2 are a root of g2, which we call beta1,
  * the other root is called gamma1
  * then beta2 and gamma2 are the roots of (Ex^2-Bx+C)
  * Note: We directly apply the transformation t1: x -> (x-alpha2)/(x-alpha1)
  */
  lc := (alpha1-alpha2)^2;
  if g22 eq 0 then
    beta2 := 1;
    if g20 eq alpha1 then
      beta1 := -g11/(2*alpha1);
      gamma1 := -2*alpha2/g11;
    else
      beta1 := -2*alpha2/g11;
      gamma1 := -g11/(2*alpha1);
    end if;
    gamma2 := (C-alpha2*B)/(C-alpha1*B);
    lc := 2*lc*alpha1*g11*(C-alpha1*B);
  else
    if g20 eq -alpha1*(alpha1-g21) then
      beta1 := -g11/(2*alpha1);
      beta2 := (g21-alpha1+alpha2)/g21;
      gamma1 := -2*alpha2/g11;
      gamma2den := B+E*(g21-2*alpha1);
      gamma2 := (B+E*(g21+g11))/gamma2den;
      lc := 2*lc*alpha1*g21*g11*gamma2den;
    else
      beta1 := -2*alpha2/g11;
      beta2 := g21/(g21+alpha1-alpha2);
      gamma1 := -g11/(2*alpha1);
      gamma2den := B + E*(g21+g11);
      gamma2 := (B+E*(g21-2*alpha2))/gamma2den;
      lc := 2*lc*alpha1*g11*(g21+alpha1-alpha2)*gamma2den;
    end if;
  end if;

  if #P eq 5 then
    /* Transformation (x-alpha2)/(x-alpha1) is applied  to the 4-torsion point P */
    a0,a1,a2,b0,b1 := Explode(P);
    mu := a1*b0-a0*b1;
    deni := 1/(alpha1*(a1 + alpha1) + a0);
    a0p :=  (alpha2*(a1 + alpha2) + a0)*deni;
    a1p :=  (a1*g11 - 2*(g10 + a0))*deni;
    aux1 := (deni*(alpha2 - alpha1))^2;
    b0p :=  aux1 * (alpha1*(alpha1*(b1*alpha2+b0) + 2*b0*alpha2) + mu*(2*alpha1+alpha2+a1) - a0*b0);
    b1p :=  aux1 * (a0*b0 - alpha1^2*(b1*alpha1 +3*b0) - (3*alpha1 +a1)*mu);
  end if;

  /* compute sqrt(beta1*beta2) as explained in Corollary 2.8 */
  beta12 :=beta1*beta2;
  if #P eq 5 then
    invsqbeta12 := (beta12*b0p^2 + lc*(a0p-beta12)*(-a1p-beta1-beta2)*a0p^2)/(b0p*(b1p*a0p-a1p*b0p)*beta12+lc*((a0p-beta12)*a0p)^2);
    //print(invsqbeta12^2 eq 1/beta12 );
  else
    invsqbeta12 := Sqrt(1/beta12);
  end if;

  //print invsqbeta12^2 eq 1/beta12;

  /* final transformation */
  alpha := invsqbeta12;
  beta := -invsqbeta12*alpha2;
  gamma := 1;
  delta := -alpha1;
  epsilon :=  (invsqbeta12*(alpha2-alpha1))^3;

  Ap := (beta1+beta2)*invsqbeta12;
  Bp := (gamma1+gamma2)*invsqbeta12;
  Cp := gamma1*gamma2*invsqbeta12^2;
  Ep := 1/(lc*alpha);

  return [alpha, beta, gamma, delta ,epsilon], [Ap,Bp,Cp,Ep];

end function;

RichelotChain := function(type2_invariants, kernel,n: P_list:=[], last_step:=false)
  /*
  * INPUT:   - type2_invariants = [A,B,C,d] defining hyperelliptic curve
  *              y^2 = (x^2-1)(x^2-A)(E*x^2-B*x+C)
  *           - kernel = [P1,P2]: Mumford coordinates of P1,P2 in defining a
  *             (2^n,2^n)-group of the Jacobian of the hyperelliptic curve
  *           - n
  *           - P_list: list of points
  *             (as in Remark 5.5)
  * OUTPUT:   - type2_invariants defining the codomain curve of the (2^n,2^n)-isogeny
  *           - Q_list: list of the image points in P_list under the isogeny
  */

  pos := 1;
  kernel_aux := kernel;
  indices := [0];
  for i := 0 to n-2 do
    //print indices;
    gap := n-1-i - indices[pos];
    if gap eq 1 then
      P4 := kernel_aux[2*pos -1];
      kernel2 := DoubleType2(kernel_aux[[2*pos-1, 2*pos]], type2_invariants);
      if not indices[pos] eq 0 or not last_step then
               Prune(~indices);
               Prune(~kernel_aux);
               Prune(~kernel_aux);
               pos -:= 1;
      end if;
    elif gap eq 2 then
      kernel4 :=  DoubleType2(kernel_aux[[2*pos-1, 2*pos]], type2_invariants);
      P4 := kernel4[1];
      kernel2 := DoubleType2(kernel4, type2_invariants);
    else
      new_ind := indices[pos] + Floor(gap/2);
      new_aux := kernel_aux[[2*pos-1,2*pos]];
      for j := 1 to Floor(gap/2) do
        new_aux := DoubleType2(new_aux, type2_invariants);
      end for;
      Append(~indices, new_ind);
      kernel_aux cat:=  new_aux;
      pos +:= 1;
      for j := 1 to Ceiling(gap/2) -1 do
        new_aux := DoubleType2(new_aux, type2_invariants);
      end for;
      P4 := new_aux[1];
      kernel2 := DoubleType2(new_aux, type2_invariants);
    end if;

    trafo, type1_invariants := FindTransformation(type2_invariants, kernel2[1], kernel2[2] : P:=P4);
    Q_lists := MumfordTransformationga1(trafo, [kernel_aux, P_list]);
    type2_invariants, Q_lists := RichelotType1(type1_invariants, Q_lists);
    kernel_aux, P_list := Explode(Q_lists);
  end for;

  if last_step then
    kernel2 := kernel_aux;
    trafo, type1_invariants := FindTransformation(type2_invariants, kernel2[1], kernel2[2]);
    Q_lists := MumfordTransformationga1(trafo, [P_list]);
    type2_invariants, P_list := RichelotType1(type1_invariants, Q_lists);
  end if;
  return type2_invariants, P_list;

end function;
