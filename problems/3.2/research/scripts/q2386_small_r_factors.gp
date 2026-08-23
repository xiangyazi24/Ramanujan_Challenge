\\ Q2386 preliminary exact candidate-prime audit.
P(n)=34*n^3+51*n^2+27*n+5;
B=vector(13);
B[1]=1; \\ b_0
B[2]=5; \\ b_1
for(n=1,11,B[n+2]=(P(n)*B[n+1]-n^3*B[n])/(n+1)^3);
print("Q2386_FACTORS_BEGIN");
for(r=0,12,b=B[r+1];print("R=",r);print("B=",b);print("DIGITS=",#Str(abs(b)));print("FACTOR=",factor(abs(b)));print("---"));
print("Q2386_FACTORS_END");
quit;
