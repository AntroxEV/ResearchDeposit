clear,clc,close all
syms r0 beta mu alpha z

ra=0.5*r0*(2+beta);
Aa=0.25*pi*r0^2*(mu*(2+alpha)^2-(2+beta)^2);
Ia=1/64*pi*r0^4*(mu^2*(2+alpha)^4-(2+beta)^4);

Runif=simplify(Aa/Ia)


A=pi*r0^2*(mu*(1+alpha*z).^2-(1+beta*z).^2);
I=0.25*pi*r0^4*(mu^2*(1+alpha*z).^4-(1+beta*z).^4);

Rtp=simplify(A/I)