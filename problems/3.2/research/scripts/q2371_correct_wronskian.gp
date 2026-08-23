\\ Q2371: correctly indexed mixed Wronskian for N_h versus R_{h-1}.
\\ Exact PARI/GP audit for gaps h=3,4,5.

default(parisizemax, 8000000000);
allocatemem(1000000000);

P(t)=34*t^3+51*t^2+27*t+5;
Q(m,x)=(2*x+2*m+1)*(3*x^2+2*x*m+m^2+3*x+m+1);
C(t)=16*t^3+24*t^2+12*t+2;

x='x;

\\ N[i+1]=N_i, R[i+1]=R_i.
N=vector(7); N[1]=0; N[2]=1;
for(m=1,5,N[m+2]=P(x+m)*N[m+1]-(x+m)^6*N[m]);
R=vector(7); R[1]=1; R[2]=Q(1,x);
for(m=2,5,R[m+1]=Q(m,x)*R[m]-(x+m)^6*R[m-1]);

\\ Actual mixed Wronskian at step m:
\\ W_m=N_{m+1}R_{m-1}-N_mR_m.
\\ delta_m=P(x+m)-Q_m=2(C(x+m)-(2(x+m)+1)lambda_x).
U=vector(6); V=vector(6);
t=x+1; U[1]=2*C(t); V[1]=2*t+1;
for(m=2,5,t=x+m;U[m]=t^6*U[m-1]+2*C(t)*N[m+1]*R[m];V[m]=t^6*V[m-1]+(2*t+1)*N[m+1]*R[m]);

audit(H) = {
  my(M=H-1,W=U[H-1]-2*x*(x+1)*V[H-1],Wdirect=N[H+1]*R[H-1]-N[H]*R[H],co,minco,pos=1,j,fac);
  print("============================================================");
  print("H=",H," M=",M);
  print("DIRECT_IDENTITY_OK=",W==Wdirect);
  print("DEG_U=",poldegree(U[M])," DEG_V=",poldegree(V[M])," DEG_W=",poldegree(W));
  print("CONTENT_U=",content(U[M])," CONTENT_V=",content(V[M])," CONTENT_W=",content(W));
  print("U=",U[M]);
  print("V=",V[M]);
  print("W=",W);
  print("FACTOR_U=",factor(U[M]));
  print("FACTOR_V=",factor(V[M]));
  fac=factor(W); print("FACTOR_W=",fac);
  print("IRRED_W=",polisirreducible(W));
  print("GCD_U_V=",gcd(U[M],V[M]));
  print("GCD_W_V=",gcd(W,V[M]));
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
  print("W_AT_MINUS2=",subst(W,x,-2));
  print("W_AT_MATE=",subst(W,x,-(H+1)/2));
  print("TWO_V_AT_MINUS2=",subst(2*V[M],x,-2));
  print("TWO_V_AT_MATE=",subst(2*V[M],x,-(H+1)/2));
  print("---")
};

print("Q2371_CORRECT_WRONSKIAN_BEGIN");
audit(3);
audit(4);
audit(5);
print("Q2371_CORRECT_WRONSKIAN_END");
quit;
