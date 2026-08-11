\\ Q7320 exact small-gap factorization probe (PARI/GP).
\\ Run: gp -q q7320_small_factor.gp
P(t)=(2*t+1)*(17*t^2+17*t+5);
H=8;
N=vector(H+1);
N[1]=0;
N[2]=1;
for(h=1,H-1,N[h+2]=P(x+h)*N[h+1]-(x+h)^6*N[h]);
print("Q7320_SMALL_FACTOR H=",H);
for(h=1,H,lhs=subst(N[h+1],x,-h-1-x);rhs=(-1)^(h-1)*N[h+1];if(lhs!=rhs,error("reflection failure at h=",h)));
for(h=2,H,if(h%2==0,L=2*x+h+1;qr=divrem(N[h+1],L);if(qr[2]!=0,error("linear factor failure at h=",h));M=qr[1];if(subst(M,x,-h-1-x)!=M,error("quotient reciprocity failure at h=",h));print("h=",h," degree(N)=",poldegree(N[h+1])," degree(M)=",poldegree(M));print("factor(N_h)=",factor(N[h+1]));print("factor(M_h)=",factor(M))));
print("Q7320_SMALL_FACTOR PASS");
quit;
