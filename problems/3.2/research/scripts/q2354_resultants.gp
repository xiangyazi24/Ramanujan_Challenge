\\ Q2354 exact fixed-gap Apéry/Racah resultants.
\\ Temporary audit script; PARI/GP, integer arithmetic only.

default(parisizemax, 8000000000);
allocatemem(1000000000);

P(t) = 34*t^3 + 51*t^2 + 27*t + 5;
Q(m,t) = (2*t + 2*m + 1) * (3*t^2 + 2*t*m + m^2 + 3*t + m + 1);

apery_gap(h,t) = {
  my(a = 1, b, c);
  if (h == 1, return(a));
  b = P(t + 1);
  if (h == 2, return(b));
  for (k = 2, h - 1,
    c = P(t + k) * b - (t + k)^6 * a;
    a = b;
    b = c;
  );
  return(b);
};

racah_gap(m,t) = {
  my(a = 1, b, c);
  if (m == 0, return(a));
  b = Q(1,t);
  if (m == 1, return(b));
  for (k = 2, m,
    c = Q(k,t) * b - (t + k)^6 * a;
    a = b;
    b = c;
  );
  return(b);
};

print("Q2354_PARIGP_BEGIN");
for (h = 3, 6,
  my(A = apery_gap(h,x), B = racah_gap(h-1,x));
  my(mate = 2*x + h + 1, Ar = A, Br = B, qrA, qrB);
  print("H=", h);
  print("DEG_N=", poldegree(A), " DEG_R=", poldegree(B));
  print("CONTENT_N=", content(A), " CONTENT_R=", content(B));
  print("REFLECT_N_OK=", subst(A,x,-x-h-1) == (-1)^(h-1)*A);
  print("REFLECT_R_OK=", subst(B,x,-x-h-1) == (-1)^(h-1)*B);
  if (h % 2 == 0,
    qrA = divrem(A,mate);
    qrB = divrem(B,mate);
    if (qrA[2] != 0 || qrB[2] != 0, error("mate division failed"));
    Ar = qrA[1]; Br = qrB[1];
    print("RAW_RESULTANT=0");
    print("COMMON_MATE_FACTOR=", mate);
  ,
    print("RAW_RESULTANT_NONZERO_EXPECTED=1");
  );
  print("DEG_N_REDUCED=", poldegree(Ar), " DEG_R_REDUCED=", poldegree(Br));
  print("CONTENT_N_REDUCED=", content(Ar), " CONTENT_R_REDUCED=", content(Br));
  my(g = gcd(Ar,Br));
  print("GCD_REDUCED=", g);
  my(res = polresultant(Ar,Br,x));
  print("REDUCED_RESULTANT_SIGN=", sign(res));
  print("REDUCED_RESULTANT_DIGITS=", #Str(abs(res)));
  print("REDUCED_RESULTANT=", res);
  my(fac = factor(abs(res)));
  print("FACTORIZATION=", fac);
  print("LARGEST_PRIME_FACTOR=", fac[matsize(fac)[1],1]);
  print("NUMBER_DISTINCT_PRIMES=", matsize(fac)[1]);
  print("---");
);
print("Q2354_PARIGP_END");
quit;
