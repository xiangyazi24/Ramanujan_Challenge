\\ Exact audit for Q2371: affine spectral split for h=3,4,5.
\\ PARI/GP, exact integer/rational arithmetic only.

default(parisizemax, 8000000000);
allocatemem(1000000000);

P(t)=34*t^3+51*t^2+27*t+5;
Q(m,x)=(2*x+2*m+1)*(3*x^2+2*x*m+m^2+3*x+m+1);
AA(t)=16*t^3+27*t^2+12*t+3;

x='x;

\\ N[i+1]=N_i, R[i+1]=R_i.
N=vector(7); N[1]=0; N[2]=1;
for(m=1,5,N[m+2]=P(x+m)*N[m+1]-(x+m)^6*N[m]);
R=vector(7); R[1]=1; R[2]=Q(1,x);
for(m=2,5,R[m+1]=Q(m,x)*R[m]-(x+m)^6*R[m-1]);

\\ U[m],V[m] use the Q2371/Q2375 normalization W_1=1.
U=vector(6); V=vector(6); U[1]=1; V[1]=0;
for(m=1,5,t=x+m+1;U[m+1]=t^6*U[m]+2*AA(t)*N[m+1]*R[m+1];V[m+1]=t^6*V[m]+(2*t-1)*N[m+1]*R[m+1]);

audit(H) = {
  my(W=U[H]-2*x*(x+1)*V[H],co,minco,pos=1,j);
  print("============================================================");
  print("H=",H);
  print("DEG_U=",poldegree(U[H])," DEG_V=",poldegree(V[H])," DEG_W=",poldegree(W));
  print("CONTENT_U=",content(U[H])," CONTENT_V=",content(V[H])," CONTENT_W=",content(W));
  print("U=",U[H]);
  print("V=",V[H]);
  print("W=",W);
  print("FACTOR_U=",factor(U[H]));
  print("FACTOR_V=",factor(V[H]));
  print("FACTOR_W=",factor(W));
  print("IRRED_W=",polisirreducible(W));
  print("GCD_U_V=",gcd(U[H],V[H]));
  print("GCD_W_V=",gcd(W,V[H]));
  print("GCD_N_R=",gcd(N[H+1],R[H]));
  print("GCD_N_W=",gcd(N[H+1],W));
  print("GCD_R_W=",gcd(R[H],W));
  print("TRIPLE_GCD=",gcd(gcd(N[H+1],R[H]),W));
  print("REAL_ROOT_COUNT_W=",polsturm(W));
  co=Vec(W); minco=co[1];
  for(j=1,#co,if(co[j]<minco,minco=co[j]);if(co[j]<=0,pos=0));
  print("ALL_W_COEFFICIENTS_POSITIVE=",pos," MIN_W_COEFFICIENT=",minco);
  print("W_AT_0=",subst(W,x,0));
  print("W_AT_MINUS1=",subst(W,x,-1));
  print("W_AT_MATE=",subst(W,x,-(H+1)/2));
  print("RATIO_DEN_AT_MATE=",subst(2*V[H],x,-(H+1)/2));
  print("---")
};

print("Q2371_AFFINE_SPLIT_BEGIN");
audit(3);
audit(4);
audit(5);
print("Q2371_AFFINE_SPLIT_END");
quit;
