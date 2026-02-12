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
m=0.1;
mu=(1+m)^2;
alphav=linspace(0.25,1,100);
betav=linspace(0,sqrt(mu),10);
% betav=[0 0.5*alpha alpha alpha*sqrt(mu)];
r0=1;
E=1;
rho=1;
AA=[]

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
        
        %
        vect=1:np*2;
        % vect([1,np*2-1])=[];
        vect([1,2,np*2-1,np*2])=[];
        Kgbc=Kg(vect,vect);
        Mgbc=Mg(vect,vect);
        
        wq=eig(Kgbc,Mgbc);
        wq=wq(~isinf(wq));
        lambda4=(4*rho*wq*z(end)^4/(E*r0^2));
        
        % lambdan=(lambda4./(1+mu)).^.25
        lambda=(lambda4).^.25;
        if ~isreal(lambda(end))
            lambda(end)=NaN;

        end
        % AA=[AA; [lambda(end) lambda(end-1) lambda(end-2)]]
        FREQ(tt,jj)=lambda(end);
    end

end


[X,Y] = meshgrid(alphav,betav);
contourf(X,Y,FREQ)
xlabel('alpha')
ylabel('beta')
