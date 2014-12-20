% Using MATLAB Nonlinear Programming Function fmincon to find the optimum solution of
% the Minimum Cost Swing Up Cart And Pendulum Problem

clear all;

dt = 0.01; % ALSO SPECIFY IN 'intfn.m', 'f.m', 'J.m'
tf = 3.6;   % ALSO SPECIFY IN 'J.m', 'intfn.m'
usat = Inf; % 25

% Initial design values and lower and upper bounds
u = 0*[0:dt:tf]; % x = [u]

ulen = length(u);

% vlb = [];
% vub = [];
vlb(1:ulen) = -usat;
vub(1:ulen) = usat;

% Optimization options
%	'MaxFunEvals', 2000 , ...
%    'TolCon', .0001, ...
options = optimset('Display','iter', ...
    'Algorithm', 'sqp', ...
	'LargeScale', 'off', ...
    'GradObj', 'on', ...
    'TolX', .000001);

% [u_opt, fval, exitflag, output, lambda] = fmincon('J', u, ...
%     [], [], [], [], vlb, vub, 'NONLCON', options);
[u_opt, fval, exitflag, output, lambda] = fmincon('J', u, ...
    [], [], [], [], vlb, vub, [], options);

u_opt
