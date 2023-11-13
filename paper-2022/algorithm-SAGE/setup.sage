from sage.misc.prandom import randrange
def RandomSymplecticGroup(B, power2):

  r1 = randint(0,2^power2-1);
  r2 = randint(0,2^power2-1);
  r3 = randint(0,2^power2-1);

  R1 = B[0] + r1*B[2] + r2*B[3];
  R2 = B[1] + r2*B[2] + r3*B[3];

  R = [R1,R2];
  R_coord = [];
  for P in R:
    [a0,a1,a2] = P[0].coefficients()
    [b0,b1] = P[1].coefficients()
    R_coord.append([a0,a1,a2,b0,b1]);
  return R_coord;
