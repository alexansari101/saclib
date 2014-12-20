function [f] = f(t,x,u)
% u is the control vector.
% x is a 2 vector [theta, theta']

% Constants
dt = 0.01;
h = 2.0;  % pend length
g = 9.81; % gravity

% index of the control to apply based on the current time.
uindx = floor(t/dt)+1;

f = [ x(2); ( g*sin(x(1)) + u(uindx)*cos(x(1)) ) / h ];

