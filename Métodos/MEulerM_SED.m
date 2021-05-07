function [t,u,v] = MEulerM_SED(f,g,a,b,n,u0,v0)
%NEULER_SED - Método de Euler para um Sistema de SED/PVI
%15/04/2021 - Tomás Silva - a2020143845@isec.pt

h = (b-a) / n;
t = a:h:b;
u = zeros(1, n+1);
v = zeros(1, n+1);

u(1) = u0;
v(1) = v0;

for i = 1:n
    u(i+1) = u(i) + (h/2) * (f(t(i), u(i), v(i)) + f(t(i+1), u(i+1), v(i+1)));
    v(i+1) = v(i) + (h/2) * (g(t(i), u(i), v(i)) + g(t(i+1), u(i+1), v(i+1)));
end



