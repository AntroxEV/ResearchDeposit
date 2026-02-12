function [Kst]=PrismaticBeam(E,I0,L)

Kst =[ (12*E*I0)/L^3,  (6*E*I0)/L^2, -(12*E*I0)/L^3,  (6*E*I0)/L^2;
  (6*E*I0)/L^2,    (4*E*I0)/L,  -(6*E*I0)/L^2,    (2*E*I0)/L;
-(12*E*I0)/L^3, -(6*E*I0)/L^2,  (12*E*I0)/L^3, -(6*E*I0)/L^2;
  (6*E*I0)/L^2,    (2*E*I0)/L,  -(6*E*I0)/L^2,    (4*E*I0)/L];
