function [x1] = intfn(u)
% Perform integration step

% Constants
dt = .01; % integration time step
tf = 3.6;
h = 2.0;  % pend length
g = 9.81; % gravity

x0 = [pi+.1, 0];

% State vector x = [ theta(t_i), theta'(t_i) ].  x0(3) = u(t_1).
fxu=@(t,x)f(t,x,u);

[T,x1] = ode45(fxu,[0:dt:tf],x0);

% Euler Integration of Dynamics: x2 = x1 + dt*( f(x1) )
% x1 = x0(1:2)' + dt*( [ x0(2); ( g*sin(x0(1)) + x0(3)*cos(x0(1)) ) / h ] );
% x1 = x1';
