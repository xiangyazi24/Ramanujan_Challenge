\\ Q2375: exact U_m,V_m split, m=1,...,5, over Z[r].

P(t)=34*t^3+51*t^2+27*t+5;
Q(m,t)=(2*t+2*m+1)*(3*t^2+2*t*m+m^2+3*t+m+1);
A(m,t)=16*(t+m)^3+27*(t+m)^2+12*(t+m)+3;
Dlt(m,t)=2*(A(m,t)-(2*(t+m)-1)*t*(t+1));

N0=0; N1=1;
N2=P(r+1)*N1-(r+1)^6*N0;
N3=P(r+2)*N2-(r+2)^6*N1;
N4=P(r+3)*N3-(r+3)^6*N2;
N5=P(r+4)*N4-(r+4)^6*N3;

R0=1; R1=Q(1,r);
R2=Q(2,r)*R1-(r+2)^6*R0;
R3=Q(3,r)*R2-(r+3)^6*R1;
R4=Q(4,r)*R3-(r+4)^6*R2;
R5=Q(5,r)*R4-(r+5)^6*R3;

U1=1; V1=0;
U2=(r+2)^6*U1+2*A(2,r)*N1*R1;
V2=(r+2)^6*V1+(2*(r+2)-1)*N1*R1;
U3=(r+3)^6*U2+2*A(3,r)*N2*R2;
V3=(r+3)^6*V2+(2*(r+3)-1)*N2*R2;
U4=(r+4)^6*U3+2*A(4,r)*N3*R3;
V4=(r+4)^6*V3+(2*(r+4)-1)*N3*R3;
U5=(r+5)^6*U4+2*A(5,r)*N4*R4;
V5=(r+5)^6*V4+(2*(r+5)-1)*N4*R4;

W1=1;
W2=(r+2)^6*W1+Dlt(2,r)*N1*R1;
W3=(r+3)^6*W2+Dlt(3,r)*N2*R2;
W4=(r+4)^6*W3+Dlt(4,r)*N3*R3;
W5=(r+5)^6*W4+Dlt(5,r)*N4*R4;

emit(name,F)={my(c=content(F),G=F/c,fac=factor(G));print("BEGIN_",name);print("DEGREE=",poldegree(F));print("CONTENT=",c);print("EXPANDED=",F);print("PRIMITIVE_FACTORIZATION=",fac);print("FACTORBACK_OK=",factorback(fac)==G);print("END_",name)};

print("Q2375_BEGIN");
print("CHECK_W1=",W1==U1-2*r*(r+1)*V1);
print("CHECK_W2=",W2==U2-2*r*(r+1)*V2);
print("CHECK_W3=",W3==U3-2*r*(r+1)*V3);
print("CHECK_W4=",W4==U4-2*r*(r+1)*V4);
print("CHECK_W5=",W5==U5-2*r*(r+1)*V5);
emit("U1",U1); emit("V1",V1);
emit("U2",U2); emit("V2",V2);
emit("U3",U3); emit("V3",V3);
emit("U4",U4); emit("V4",V4);
emit("U5",U5); emit("V5",V5);
print("Q2375_END");
quit;
