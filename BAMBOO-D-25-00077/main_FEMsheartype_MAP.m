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
gamma=0
m=0.1;

mu=(1+m)^2;

alphav=linspace(0.1,1,100);
betav=linspace(0,0.5,100);
r0=1;
E=1;
rho=1;

for jj = 1:length(alphav)
    alpha=alphav(jj);
    for tt = 1:length(betav)
        beta=betav(tt);

        A=pi*r0^2*(mu*(1+alpha*z).^2-(1+beta*z).^2);
        I=0.25*pi*r0^4*(mu^2*(1+alpha*z).^4-(1+beta*z).^4);
    
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
    vect=1:np*2;
    vect([2,np*2-1,np*2])=[];
    
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
    if ~isreal(lambda(1))
        lambda(1)=NaN;

    end

    AA=[lambda(1) lambda(2) lambda(3)];


    
    % tau=zeros(length(vect),1);
    % indt=diag(Mgbc)~=0;
    % tau(indt)=1;
    % 
    % for hh = 1 : length(indw)
    % gammam(hh) = PHIm(:,hh).'*Mgbc*tau;
    % meff(hh)=gammam(hh)^2/(PHIm(:,hh).'*Mgbc*PHIm(:,hh));
    % end
    
    % Mtot=sum(diag(Mgbc));
    % Meff=meff./Mtot
    % sum(meff)/Mtot
    
    
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
    if ~isreal(lambda(1))
        lambda(1)=NaN;

    end
    BB=[lambda(1) lambda(2) lambda(3)];
    
    % tau=zeros(length(vect),1);
    % indt=diag(Mgbc2)~=0;
    % tau(indt)=1;
    % 
    % for hh = 1 : length(indw)
    % gamma(hh) = PHIm(:,hh).'*Mgbc2*tau;
    % meff2(hh)=gamma(hh)^2/(PHIm(:,hh).'*Mgbc2*PHIm(:,hh));
    % end
    % 
    % Mtot2=sum(diag(Mgbc2));
    % Meff=meff2./Mtot2
    % sum(meff2)/Mtot2
    % 
    % f=zeros(length(vect),1);
    % f(1)=1E10;
    % u2=Kgbc2\f;
    % u1=Kgbc\f;
    
    
    % 1-u1(1)/u2(1);
    
    err=1-AA(1)./BB(1);
    
    % errm=[ 1-meff(1)/meff2(1) 1-meff(2)/meff2(2) 1-meff(3)/meff2(3)];

ERR(tt,jj)=abs(err);
    end
end

[X,Y] = meshgrid(alphav,betav);
contourf(X,Y,ERR,linspace(0,0.3,16),"ShowText",true,"FaceAlpha",0.75)
xlabel('alpha')
ylabel('beta')