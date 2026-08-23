\\ Q2354 exact fixed-gap Apéry/Racah resultants.
\\ Temporary audit script; PARI/GP, integer arithmetic only.

default(parisizemax, 8000000000);
allocatemem(1000000000);

P(t) = 34*t^3 + 51*t^2 + 27*t + 5;
Q(m,t) = (2*t + 2*m + 1) * (3*t^2 + 2*t*m + m^2 + 3*t + m + 1);

apery_gap(H,T) = {
  my(a = 1, b, c);
  if (H == 1, return(a));
  b = P(T + 1);
  if (H == 2, return(b));
  for(k = 2, H - 1, c = P(T + k) * b - (T + k)^6 * a; a = b; b = c);
  b
};

racah_gap(M,T) = {
  my(a = 1, b, c);
  if (M == 0, return(a));
  b = Q(1,T);
  if (M == 1, return(b));
  for(k = 2, M, c = Q(k,T) * b - (T + k)^6 * a; a = b; b = c);
  b
};

audit(H) = {
  my(A = apery_gap(H,X), B = racah_gap(H-1,X));
  my(mate = 2*X + H + 1, Ar = A, Br = B, qrA, qrB, g, res, fac, nf);
  print("H=", H);
  print("DEG_N=", poldegree(A), " DEG_R=", poldegree(B));
  print("CONTENT_N=", content(A), " CONTENT_R=", content(B));
  print("REFLECT_N_OK=", subst(A,X,-X-H-1) == (-1)^(H-1)*A);
  print("REFLECT_R_OK=", subst(B,X,-X-H-1) == (-1)^(H-1)*B);
  if(H % 2 == 0, qrA = divrem(A,mate); qrB = divrem(B,mate); if(qrA[2] != 0 || qrB[2] != 0, error("mate division failed")); Ar = qrA[1]; Br = qrB[1]; print("RAW_RESULTANT=0"); print("COMMON_MATE_FACTOR=", mate), print("RAW_RESULTANT_NONZERO_EXPECTED=1"));
  print("DEG_N_REDUCED=", poldegree(Ar), " DEG_R_REDUCED=", poldegree(Br));
  print("CONTENT_N_REDUCED=", content(Ar), " CONTENT_R_REDUCED=", content(Br));
  g = gcd(Ar,Br);
  print("GCD_REDUCED=", g);
  res = polresultant(Ar,Br,X);
  print("REDUCED_RESULTANT_SIGN=", sign(res));
  print("REDUCED_RESULTANT_DIGITS=", #Str(abs(res)));
  print("REDUCED_RESULTANT=", res);
  fac = factor(abs(res));
  nf = matsize(fac)[1];
  print("FACTORIZATION=", fac);
  print("LARGEST_PRIME_FACTOR=", fac[nf,1]);
  print("NUMBER_DISTINCT_PRIMES=", nf);
  print("---")
};

print("Q2354_PARIGP_BEGIN");
audit(3);
audit(4);
audit(5);
audit(6);
print("Q2354_PARIGP_END");
quit;
