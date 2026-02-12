clear,clc,close all
syms('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','c13','c14','c15','c16','wi','z')
syms('E','L','A','I','rho')
syms('alpha','lambda')
%% constant
% lambda1=L*(rho*A*wi^2/(E*I))^(1/4);
k1(z)=lambda*z;
vmp(z)=c1*sin(k1)+c2*cos(k1)+c3*sinh(k1)+c4*cosh(k1);

Dvmp = diff(vmp,z,1)./L;
DDvmp=diff(vmp,z,2)./L^2;

bc1=vmp(0)==0; % v1
bc2=Dvmp(0)==0; %,phi1,...
bc3=vmp(1)==0; % v2
bc4=Dvmp(1)==0; %,phi12,...
SEG1=[bc1,bc2,bc3,bc4];
[B1,b] = equationsToMatrix(SEG1,c1, c2, c3, c4);

M1(z) =  DDvmp;
V1(z)=diff(M1,z,1)./L;

sc1=V1(0) == 0; %applied V1.
sc2=-M1(0) == 0; %applied M1
sc3=-V1(1) == 0; %applied V1.
sc4=M1(1) == 0; %applied M1
SSEG1=[sc1,sc2,sc3,sc4];
[C1,b2] = equationsToMatrix(SSEG1,c1, c2, c3, c4);

Kdyn1=simplify(C1/B1) %(Cn Bn^-1)  
chr = latex(Kdyn1);

Kdyn1=E*I *Kdyn1