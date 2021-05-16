function [t,u,v] = NRK4_SED(f,g,a,b,n,u0,v0)

%RK4 - Método de Runge-Kutta 4 para resolução numérica de SED
%   u'=f(t,u,v), t=[a,b], u(a)=u0
%   v'=f(t,u,v), t=[a,b], v(a)=v0
%   u(i+1) = u(i)+(k1u+2*k2u+2*k3u+k4u)/6;
%   v(i+1) = v(i)+(k1v+2*k2v+2*k3v+k4v)/6;

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

for i = 1:n % Aplicar o método de RK4 (iteração)
    k1u = h*f(t(i), u(i), v(i)); % Cálculo do valor de k1 para a função u
    k1v = h*g(t(i), u(i), v(i)); % Cálculo do valor de k1 para a função v
    
    k2u = h*f(t(i) + h/2, u(i) + k1u/2, v(i) + k1v/2); % Cálculo do valor de k2 para a função u
    k2v = h*g(t(i) + h/2, u(i) + k1u/2, v(i) + k1v/2); % Cálculo do valor de k2 para a função v
    
    k3u = h*f(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2); % Cálculo do valor de k3 para a função u
    k3v = h*g(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2); % Cálculo do valor de k3 para a função v
    
    k4u = h*f(t(i+1), u(i) + k3u, v(i) + k3v); % Cálculo do valor de k4 para a função u
    k4v = h*g(t(i+1), u(i) + k3u, v(i) + k3v); % Cálculo do valor de k4 para a função v
    
    u(i+1) = u(i) + (k1u + 2*k2u + 2*k3u + k4u) / 6; % Cálculo para a função u
    v(i+1) = v(i) + (k1v + 2*k2v + 2*k3v + k4v) / 6; % Cálculo para a função v
end
end

