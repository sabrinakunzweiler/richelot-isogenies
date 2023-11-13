function FromJacToJac(h, kernel, points, a)
  D11,D12,D21,D22 := Explode(kernel);

  R<x> := Parent(h);
  Fp2 := BaseRing(R);

  J := Jacobian(HyperellipticCurve(h));
  D1 := elt<J | D11, D12>;
  D2 := elt<J | D21, D22>;

  G1 := (2^(a-1)*D1)[1];
  G2 := (2^(a-1)*D2)[1];
  G3 := h div (G1*G2);

  delta := Matrix(Fp2, 3, 3, [Coefficient(G1, 0), Coefficient(G1, 1), Coefficient(G1, 2),
                              Coefficient(G2, 0), Coefficient(G2, 1), Coefficient(G2, 2),
                              Coefficient(G3, 0), Coefficient(G3, 1), Coefficient(G3, 2)]);
  delta := Determinant(delta)^(-1);

  H1 := delta*(Derivative(G2)*G3 - G2*Derivative(G3));
  H2 := delta*(Derivative(G3)*G1 - G3*Derivative(G1));
  H3 := delta*(Derivative(G1)*G2 - G1*Derivative(G2));

  hnew := H1*H2*H3;
  Jnew := Jacobian(HyperellipticCurve(hnew));

  // now compute image points
  // first the point [D11, D12]:

  u0 := Coefficient(D11, 0);
  u1 := Coefficient(D11, 1);
  v0 := Coefficient(D12, 0);
  v1 := Coefficient(D12, 1);
  S<x1,y1,y2,x2> := PolynomialRing(Fp2, 4);
  pr := hom<S -> R | 0, 0, 0, x>;

  eq1 := x1^2 + u1*x1 + u0;
  eq2 := v1*x1 + v0 - y1;
  eq3 := Evaluate(G1, x1)*Evaluate(H1, x2) + Evaluate(G2, x1)*Evaluate(H2, x2);
  eq4 := y1*y2 - Evaluate(G1, x1)*Evaluate(H1, x2)*(x1 - x2);
  eq5 := y1^2 - Evaluate(h, x1);
  I := Ideal([eq1, eq2, eq3, eq4, eq5]);
  G := GroebnerBasis(I); // last two are in non-reduced Mumford form: y2 + cubic(x2), quartic(x2)
  unew := pr(G[#G]);
  vnew := -pr(G[#G-1]);
  // sanity check: (vnew^2 - hnew) mod unew;
  imD1 := elt<Jnew | pr(G[#G]), -pr(G[#G-1])>;

  // now same for the point [D21, D22]:

  u0 := Coefficient(D21, 0);
  u1 := Coefficient(D21, 1);
  v0 := Coefficient(D22, 0);
  v1 := Coefficient(D22, 1);
  eq1 := x1^2 + u1*x1 + u0;
  eq2 := v1*x1 + v0 - y1;
  I := Ideal([eq1, eq2, eq3, eq4, eq5]);
  G := GroebnerBasis(I);
  unew := pr(G[#G]);
  vnew := -pr(G[#G-1]);
  imD2 := elt<Jnew | pr(G[#G]), -pr(G[#G-1])>;
  kernel_new := [imD1[1], imD1[2], imD2[1], imD2[2]];

  points_new := [];
  for P in points do
    P1,P2 := Explode(P);
    u0 := Coefficient(P1, 0);
    u1 := Coefficient(P1, 1);
    v0 := Coefficient(P2, 0);
    v1 := Coefficient(P2, 1);
    eq1 := x1^2 + u1*x1 + u0;
    eq2 := v1*x1 + v0 - y1;
    I := Ideal([eq1, eq2, eq3, eq4, eq5]);
    G := GroebnerBasis(I);
    unew := pr(G[#G]);
    vnew := -pr(G[#G-1]);
    imP := elt<Jnew | pr(G[#G]), -pr(G[#G-1])>;
    points_new := Append(points_new, [imP[1],imP[2]]);
  end for;

  return hnew, kernel_new, points_new;
end function;

function IsogenyChainCD(kernel,points,h,a)
  for i in [0..a-1] do
    h, kernel, points := FromJacToJac(h, kernel, points, a-i);
    // print "order 2^", a - i - 1, "on hyp curve", h;
  end for;
  return h, points;
end function;
