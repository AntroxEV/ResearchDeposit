clear,clc,close all
% FEM formulation for checking tapered beams
np=25;
z=linspace(0,1,np);
dz=z(2)-z(1);
% r0=1;
% s0=0.01;
% c0=0.1;
% c1=0.2;
% A=pi*((r0+s0+c0*x).^2-(r0+ci*x).^2);
% I=0.25*pi*((r0+s0+c0*x).^4-(r0+ci*x).^4);
%% BCs
vect=1:np*2;
% vect([1,np*2-1])=[]; %pinned pinned
vect([1,2,np*2-1,np*2])=[]; %fixed fixed
% vect([np*2-1,np*2])=[]; %fixed - free
% vect([2,np*2-1,np*2])=[]; %fixed - slider
%%
gamma=0
m=0.1; %0.1 0.25 0.4
alphav=[0.25; 0.5; 0.75;1];
alphav=[1];
mu=(1+m)^2;
r0=1;
E=1;
rho=1;
AA=[];
BB=[];
CC=[];
for jj=1:length(alphav)
    alpha=alphav(jj);
    betav=[0 0.25*alpha 0.5*alpha  0.75*alpha alpha alpha*sqrt(mu)];
    for tt=1:length(betav)
    beta=betav(tt);
    A=pi*r0^2*(mu*(1+alpha*z).^2-(1+beta*z).^2);
    I=0.25*pi*r0^4*(mu^2*(1+alpha*z).^4-(1+beta*z).^4);
%% FEM static + lumped mass matrix
    Kg=zeros(np*2,np*2);
    Mg=zeros(np*2,np*2);
        
    for hh=1:np-1
        Ii(hh)=(I(hh)+I(hh+1))/2;
        Ai(hh)=(A(hh)+A(hh+1))/2;
        [Kst]=PrismaticBeam(E,Ii(hh),dz);
        ind=(hh-1)*2+1:2*hh+2;
        Kg(ind,ind)=Kg(ind,ind)+Kst;
        Mi=zeros(4,4);
        Mi(1,1)=rho*Ai(hh)*dz/2;
        Mi(3,3)=rho*Ai(hh)*dz/2;
        Mg(ind,ind)=Mg(ind,ind)+Mi;
    end
    % Mg(end-1,end-1)=Mg(end-1,end-1)+gamma*pi;
    Mg(1,1)=Mg(1,1)+gamma*pi;
    
    Kgbc=Kg(vect,vect);
    Mgbc=Mg(vect,vect);
    
    % indt=diag(Mgbc)~=0;
    % Kgbc=Kgbc(indt,indt);
    % Mgbc=Mgbc(indt,indt);
    
    [PHI,WQ]=eig(Kgbc,Mgbc);
    wq=diag(WQ);
    indinf=isinf(wq);
    wq=wq(~isinf(wq));
    PHI(:,indinf)=[];
    [wqs,indw] = sort(wq);
    PHIs=PHI(:,indw);
    for hh =1 : length(indw)
    PHIm(:,hh)=PHIs(:,hh)./sqrt(PHIs(:,hh).'*Mgbc*PHIs(:,hh));
    end
    
    lambda4=(4*rho*wqs*z(end)^4/(E*r0^2));
    
    % lambdan=(lambda4./(1+mu)).^.25
    lambda=(lambda4).^.25;
    
    AA=[AA;[lambda(1) lambda(2) lambda(3)]];
    
tau=zeros(length(vect),1);
indt=diag(Mgbc)~=0;
tau(indt)=1;

for zz = 1 : length(indw)
Lgamma(zz) = PHIm(:,zz).'*Mgbc*tau;
meff(zz)=Lgamma(zz)^2/(PHIm(:,zz).'*Mgbc*PHIm(:,zz));
end

Mtot(tt,jj)=sum(diag(Mgbc));
Meff(tt,:,jj)=meff./Mtot(tt,jj);
Meffcheck(tt,jj)=sum(meff)/Mtot(tt,jj);
    
    
    %% UNIFORM 
    Kg=zeros(np*2,np*2);
    Mg=zeros(np*2,np*2);
        
    % Iim=(I(1)+I(end))/2;
    % Aim=(A(1)+A(end))/2;
    
    
    ra=0.5*r0*(2+beta);
    Aa=0.25*pi*r0^2*(mu*(2+alpha)^2-(2+beta)^2);
    Ia=1/64*pi*r0^4*(mu^2*(2+alpha)^4-(2+beta)^4);
    
    for hh=1:np-1
        [Kst]=PrismaticBeam(E,Ia,dz);
        ind=(hh-1)*2+1:2*hh+2;
        Kg(ind,ind)=Kg(ind,ind)+Kst;
        Mi=zeros(4,4);
        Mi(1,1)=rho*Aa*dz/2;
        Mi(3,3)=rho*Aa*dz/2;
        Mg(ind,ind)=Mg(ind,ind)+Mi;
    end
    %
    % Mg(end-1,end-1)=Mg(end-1,end-1)+gamma*pi;
    Mg(1,1)=Mg(1,1)+gamma*pi;
    
    Kgbc2=Kg(vect,vect);
    Mgbc2=Mg(vect,vect);
    
    
    [PHI,WQ]=eig(Kgbc2,Mgbc2);
    wq=diag(WQ);
    indinf=isinf(wq);
    wq=wq(~isinf(wq));
    PHI(:,indinf)=[];
    [wqs,indw] = sort(wq);
    PHIs=PHI(:,indw);
    for hh =1 : length(indw)
    PHIm(:,hh)=PHIs(:,hh)./sqrt(PHIs(:,hh).'*Mgbc2*PHIs(:,hh));
    end
    
    lambda4=(4*rho*wqs*z(end)^4/(E*r0^2));
    
    % lambdan=(lambda4./(1+mu)).^.25
    lambda=(lambda4).^.25;
    
    BB=[BB;[lambda(1) lambda(2) lambda(3)]];
    
    tau=zeros(length(vect),1);
    indt=diag(Mgbc2)~=0;
    tau(indt)=1;
    
    for zz = 1 : length(indw)
    Lgamma(zz) = PHIm(:,zz).'*Mgbc2*tau;
    meff2(zz)=Lgamma(zz)^2/(PHIm(:,zz).'*Mgbc2*PHIm(:,zz));
    end
    
    Mtot2(tt,jj)=sum(diag(Mgbc2));
    Meff2(tt,:,jj)=meff2./Mtot2(tt,jj);
    Meffcheck(tt,jj)=sum(meff)/Mtot2(tt,jj);
    % 
    f=zeros(length(vect),1);
    f(1)=1;
    u2=Kgbc2\f;
    k2(tt,jj)=f(1)/u2(1);
    u1=Kgbc\f;
    k1(tt,jj)=f(1)/u1(1);
    
%% FEM dynamic
for xx=1:3
lq=AA(end,xx)*0.9:0.001:AA(end,xx)*1.1;
fx=0;
for nn = 1 :length(lq)
    Kg=zeros(np*2,np*2);
    for hh=1:np-1
        Ii(hh)=(I(hh)+I(hh+1))/2;
        Ai(hh)=(A(hh)+A(hh+1))/2;
        [Kdyn] = DynamicConstBeam(E,dz,Ai(hh),Ii(hh),rho,lq(nn).^2/2);
        ind=(hh-1)*2+1:2*hh+2;
        Kg(ind,ind)=Kg(ind,ind)+Kdyn;
    end
    Kgbc=Kg(vect,vect);
    fx(nn)=det(Kgbc);
end
wi1=lq(abs(fx)==min(abs(fx))).^2/2;

lambda4=(4*rho*wi1.^2*z(end)^4/(E*r0^2));
    
    % lambdan=(lambda4./(1+mu)).^.25
lambdai(xx)=(lambda4).^.25;
end
    
    CC=[CC;[lambdai(1) lambdai(2) lambdai(3)]];
    % 1-u1(1)/u2(1)
    
    err=1-BB.^2./AA.^2;
    
    % errm=[ 1-meff(1)/meff2(1) 1-meff(2)/meff2(2) 1-meff(3)/meff2(3)]
    end

end
err
AA
BB
