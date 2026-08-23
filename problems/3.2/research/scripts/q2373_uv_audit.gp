\\ Q2373: exact U_h,V_h split, polynomial gcds, and finite-field audits.

P(t)=34*t^3+51*t^2+27*t+5;
Q(m,t)=(2*t+2*m+1)*(3*t^2+2*t*m+m^2+3*t+m+1);
Au(t)=16*t^3+27*t^2+12*t+3;

N0=0; N1=1;
N2=P(r+1)*N1-(r+1)^6*N0;
N3=P(r+2)*N2-(r+2)^6*N1;
N4=P(r+3)*N3-(r+3)^6*N2;
N5=P(r+4)*N4-(r+4)^6*N3;

R0=1; R1=Q(1,r);
R2=Q(2,r)*R1-(r+2)^6*R0;
R3=Q(3,r)*R2-(r+3)^6*R1;
R4=Q(4,r)*R3-(r+4)^6*R2;

U1=1; V1=0;
U2=(r+2)^6*U1+2*Au(r+2)*N1*R1;
V2=(r+2)^6*V1+(2*(r+2)-1)*N1*R1;
U3=(r+3)^6*U2+2*Au(r+3)*N2*R2;
V3=(r+3)^6*V2+(2*(r+3)-1)*N2*R2;
U4=(r+4)^6*U3+2*Au(r+4)*N3*R3;
V4=(r+4)^6*V3+(2*(r+4)-1)*N3*R3;
U5=(r+5)^6*U4+2*Au(r+5)*N4*R4;
V5=(r+5)^6*V4+(2*(r+5)-1)*N4*R4;

W1=1;
W2=(r+2)^6*W1+2*(Au(r+2)-(2*(r+2)-1)*r*(r+1))*N1*R1;
W3=(r+3)^6*W2+2*(Au(r+3)-(2*(r+3)-1)*r*(r+1))*N2*R2;
W4=(r+4)^6*W3+2*(Au(r+4)-(2*(r+4)-1)*r*(r+1))*N3*R3;
W5=(r+5)^6*W4+2*(Au(r+5)-(2*(r+5)-1)*r*(r+1))*N4*R4;

NN(h)=if(h==3,N3,if(h==4,N4,N5));
RR(h)=if(h==3,R2,if(h==4,R3,R4));
UU(h)=if(h==3,U3,if(h==4,U4,U5));
VV(h)=if(h==3,V3,if(h==4,V4,V5));
WW(h)=if(h==3,W3,if(h==4,W4,W5));
ev(F,a,p)=lift(subst(F,r,Mod(a,p)));

emitpoly(name,F)={
  my(c=content(F),G=F/c,fac=factor(G));
  print("BEGIN_",name);
  print("DEGREE=",poldegree(F));
  print("CONTENT=",c);
  print("EXPANDED=",F);
  print("PRIMITIVE_FACTORIZATION=",fac);
  print("FACTORBACK_OK=",factorback(fac)==G);
  print("END_",name)
};

gcdres(h)={
  my(gg=gcd(NN(h),VV(h)));
  print("h=",h," GCD_N_V=",gg);
  if(h==3,
    my(res=polresultant(NN(h),VV(h)));
    print("h=3 RESULTANT=",res);
    print("h=3 RESULTANT_FACTORIZATION=",factor(res))
  )
};

auditHP(h,p)={
  my(F=NN(h),G=VV(h),H=RR(h),J=UU(h));
  for(a=0,p-1,
    my(nf=ev(F,a,p),vg=ev(G,a,p));
    if(nf==0 && vg==0,
      nvcount++;
      if(a+h<p,nvrange++);
      print("NV_COMMON p=",p," h=",h," r=",a," nonwrap=",a+h<p," mate=",(2*a+h+1)%p==0)
    );
    if(nf==0 && ev(H,a,p)==0,
      my(wz=ev(WW(h),a,p),eq=lift(Mod(ev(J,a,p)-2*a*(a+1)*vg,p)),mate=((2*a+h+1)%p==0));
      commoncount++;
      if(vg==0,commonvzero++);
      if(wz==0,nrwzero++,nrwnonzero++);
      if(h==3 || vg==0 || wz==0 || !mate,
        print("NR_COMMON p=",p," h=",h," r=",a," V=",vg," U=",ev(J,a,p)," W=",wz," splitResidual=",eq," nonwrap=",a+h<p," mate=",mate)
      )
    )
  )
};

auditP(p)={for(h=3,5,auditHP(h,p))};

print("Q2373_EXACT_BEGIN");
print("SPLIT_CHECK_1=",W1==U1-2*r*(r+1)*V1);
print("SPLIT_CHECK_2=",W2==U2-2*r*(r+1)*V2);
print("SPLIT_CHECK_3=",W3==U3-2*r*(r+1)*V3);
print("SPLIT_CHECK_4=",W4==U4-2*r*(r+1)*V4);
print("SPLIT_CHECK_5=",W5==U5-2*r*(r+1)*V5);
emitpoly("V3",V3);
emitpoly("V4",V4);
emitpoly("V5",V5);

print("GCD_AND_RESULTANT_AUDIT");
for(h=3,5,gcdres(h));

print("FINITE_FIELD_SCAN_BEGIN");
nvcount=0; nvrange=0; commoncount=0; commonvzero=0; nrwzero=0; nrwnonzero=0;
forprime(p=3,500,auditP(p));
print("NV_TOTAL=",nvcount);
print("NV_NONWRAP_TOTAL=",nvrange);
print("NR_TOTAL=",commoncount);
print("NR_WITH_V_ZERO=",commonvzero);
print("NR_WITH_W_ZERO=",nrwzero);
print("NR_WITH_W_NONZERO=",nrwnonzero);
print("FINITE_FIELD_SCAN_END");
print("Q2373_EXACT_END");
quit;