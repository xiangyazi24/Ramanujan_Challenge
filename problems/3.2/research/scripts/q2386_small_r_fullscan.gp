\\ Q2386 exhaustive folded-half scan for every admissible p | b_r, 1 <= r <= 8.
P(n)=34*n^3+51*n^2+27*n+5;
Bint=vector(9);Bint[1]=1;Bint[2]=5;for(n=1,7,Bint[n+2]=(P(n)*Bint[n+1]-n^3*Bint[n])/(n+1)^3);
audit(R,PP)={my(M=(PP-1)/2,BZ=List(),KZ=List(),TZ=List(),bm1=1%PP,b0=5%PP,km1=1%PP,k0=(1+2*R*(R+1))%PP,bn,kn,idx,C);if(bm1==0,listput(BZ,0));if(km1==0,listput(KZ,0));if(M>=1,if(b0==0,listput(BZ,1));if(k0==0,listput(KZ,1));if(1>R&&b0==0&&k0==0,listput(TZ,1)));for(n=1,M-1,idx=n+1;bn=lift(Mod(P(n)*b0-n^3*bm1,PP)/Mod((n+1)^3,PP));C=(2*n+1)*(n*(n+1)+1+2*R*(R+1));kn=lift(Mod(C*k0-n^3*km1,PP)/Mod((n+1)^3,PP));if(bn==0,listput(BZ,idx));if(kn==0,listput(KZ,idx));if(idx>R&&bn==0&&kn==0,listput(TZ,idx));bm1=b0;b0=bn;km1=k0;k0=kn);print("R=",R," P=",PP," M=",M," B_R_MOD=",Bint[R+1]%PP);print("B_ZERO_FOLDED=",Vec(BZ));print("K_ZERO_FOLDED=",Vec(KZ));print("NONMATE_TRIPLE=",Vec(TZ));print("TRIPLE_COUNT=",#TZ);print("---")};
auditR(R)={my(F=factor(abs(Bint[R+1])));for(j=1,matsize(F)[1],my(PP=F[j,1]);if(PP>=2*R+1,audit(R,PP)))};
print("Q2386_FULLSCAN_BEGIN");
for(R=1,8,auditR(R));
print("Q2386_FULLSCAN_END");
quit;
