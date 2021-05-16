function [t,u,v] = NRK2_SED(f,g,a,b,n,u0,v0)

%RK2 - Método de Runge-Kutta 2 para resolução numérica de SED
%   u'=f(t,u,v), t=[a,b], u(a)=u0
%   v'=g(t,u,v), t=[a,b], v(a)=v0
%   u(i+1)=u(i)+ 1/2*(k1u+k2u)
%   v(i+1)=v(i)+ 1/2*(k1v+k2v)

%INPUT:
%   f - primeira equação diferencial
%   g - segunda equação diferencial
%   [a,b] - intervalo de valores da variável independente t
%   n - núnmero de subintervalos ou iterações do método
%   u0 - aproximação inicial u(a)=u0
%   v0 - aproximação inicial v(a)=v0

%OUTPUT:
%   t - vetor com os valores que a varíavel independente pode tomar
%   u - vetor das soluções aproximadas do SED para a função u
%   v - vetor das soluções aproximadas do SED para a função v

%   15/05/2021  Tomás Silva  a2020143845@isec.pt
%   15/05/2021  Tomás Pinto  a2020144067@isec.pt
%   15/05/2021  Francisco Mendes  a2020143982@isec.pt

h = (b-a)/n; % Amplitude de cada subintervalo
t = a:h:b; % Criar vetor que vai de "a" a "b" com step de "h"
u = zeros(1,n+1); % Alocação de memória para a função u
v = zeros(1,n+1); % Alocação de memória para a função v
u(1) = u0; % Definir o primeiro elemento do array para a função u
v(1) = v0; % Definir o primeiro elemento do array para a função v

for i = 1:n % Aplicar o método de RK2 (iteração)
    
    k1u = h*f(t(i), u(i), v(i)); % Cálculo do valor de k1 para a função u
    k1v = h*g(t(i), u(i), v(i)); % Cálculo do valor de k1 para a função v
    
    k2u = h*f(t(i+1), u(i) + k1u, v(i) + k1v); % Cálculo do valor de k2 para a função u
    k2v = h*g(t(i+1), u(i) + k1u, v(i) + k1v); % Cálculo do valor de k2 para a função v
    
    u(i+1) = u(i) + (k1u + k2u) / 2; % Cálculo para a função u
    v(i+1) = v(i) + (k1v + k2v) / 2; % Cálculo para a função v
end
end

