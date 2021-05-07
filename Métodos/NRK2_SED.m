function [t,u,v] = NRK2_SED(f,g,a,b,n,u0,v0)
%NEULERSED Método de RK2 para um Sistema de SED/PVI
%   Detailed explanation goes here
%   15/04/2021  Arménio Correia   armenioc@isec.pt

h = (b-a)/n;
t = a:h:b;
u = zeros(1,n+1);
v = zeros(1,n+1);
u(1) = u0;
v(1) = v0;
for i = 1:n
    k1u = h*f(t(i), u(i), v(i));
    k1v = h*g(t(i), u(i), v(i));
    
    k2u = h*f(t(i+1), u(i) + k1u, v(i) + k1v);
    k2v = h*g(t(i+1), u(i) + k1u, v(i) + k1v);
    
    u(i+1) = u(i) + (k1u + k2u) / 2;
    v(i+1) = v(i) + (k1v + k2v) / 2;   
end
end

