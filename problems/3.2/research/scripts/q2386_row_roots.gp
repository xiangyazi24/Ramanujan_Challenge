\\ Q2386 exact fixed-r row-root audit over candidate primes p | b_r.
P(n)=34*n^3+51*n^2+27*n+5;
rowpoly(R,T)={my(f=factorial(R),q=1,A=0,c);for(k=0,R,if(k>0,q=q*(T-(k-1)*k));c=f\factorial(k);A=A+binomial(R,k)*binomial(R+k,k)*c^2*q);A};
B=vector(13);B[1]=1;B[2]=5;for(n=1,11,B[n+2]=(P(n)*B[n+1]-n^3*B[n])\((n+1)^3));
audit(R,PP)={my(A=rowpoly(R,x),FF=factormod(A,PP),L=List(),roots,M=(PP-1)\2,inv2=(PP+1)\2,extra=0);for(ii=1,matsize(FF)[1],my(ff=FF[ii,1]);if(poldegree(ff)==1,listput(L,lift(-polcoeff(ff,0)/polcoeff(ff,1)))));roots=Vec(L);print("R=",R," P=",PP," B_R_MOD=",B[R+1]%PP," DEG_A=",poldegree(A)," LINEAR_ROOTS=",#roots," FACTOR_DEGREES=",vector(matsize(FF)[1],ii,poldegree(FF[ii,1])));for(i=1,#roots,my(t=roots[i]%PP,dlt=(1+4*t)%PP,leg=kronecker(dlt,PP),s1,s2,sf);if(dlt==0,sf=M;if(sf!=R,extra++);print("ROOT_T=",t," GRID_S=",sf," SELF=",sf==R),if(leg==1,my(d=lift(sqrt(Mod(dlt,PP))));s1=((d-1)*inv2)%PP;s2=PP-1-s1;sf=min(s1,s2);if(sf!=R,extra++);print("ROOT_T=",t," GRID_S=",sf," SELF=",sf==R),print("ROOT_T=",t," NONGRID=1"))));print("EXTRA_GRID_ROOTS=",extra);print("---")};
auditR(R)={my(F=factor(abs(B[R+1])));for(j=1,matsize(F)[1],my(PP=F[j,1]);if(PP>=2*R+1,audit(R,PP)))};
print("Q2386_ROOTS_BEGIN");
for(R=1,12,auditR(R));
print("Q2386_ROOTS_END");
quit;
