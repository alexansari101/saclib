function [x1] = intfn(x0,u)
% Perform integration from vector of controls, u

global dt tf

% State vector x = [ theta(t_i), theta'(t_i) ].  x0(3) = u(t_1).
fxu=@(t,x)f(t,x,u);

[T,x1] = ode45(fxu,[0:dt:tf],x0);

% x cannot contain x0 u0
x1 = x1(2:end,:);