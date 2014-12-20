% Using MATLAB Nonlinear Programming Function fmincon to find the optimum solution of
% the Minimum Cost Swing Up Cart And Pendulum Problem

clear all;

global dt tf xlen h x0 derivC

tic
dt = 0.003;  
tf = 6;   
h = 2.0;    % pend length
g = 9.81;   % gravity
usat = 25;  % INF

x0 = [pi, 0, 0];                % x0 = [theta theta' u0]
u0 = 0*[0:dt:tf];               % initial u
xtrue = intfn(x0(1:2),u0);      % x = [theta theta']
xlen = length(xtrue(:,1))*3;
x(1:3:xlen) = xtrue(:,1)';
x(2:3:xlen) = xtrue(:,2)';
x(3:3:xlen) = u0(2:end);

% vlb = [];
% vub = [];
%
vlb(1:3:xlen) = -Inf;
vlb(2:3:xlen) = -Inf;
vlb(3:3:xlen) = -usat;
vub(1:3:xlen) = Inf;
vub(2:3:xlen) = Inf;
vub(3:3:xlen) = usat;
%
derivC = zeros(2*xlen/3, xlen);
derivC(1,1) = -1;
derivC(2,2) = -1;
dc1 = [1 dt 0 -1]; % rows [3 : 1+xlen/3]
for i = [3:1+xlen/3]
  derivC(i,3*(i-3)+1:3*(i-3)+4) = dc1;
end
rb = 2+xlen/3;
re = 2*xlen/3;
for i = [rb:re]
  derivC(i,3*(i-rb)+1:3*(i-rb)+5) = ... 
      [0 1 0 0 -1];
end
%
options = optimset('Display','iter', ...
    'Algorithm', 'sqp', ...
	'LargeScale', 'off', ...
    'GradObj', 'on', ...
    'GradConstr','off', ...
    'MaxIter', 50, ...
    'MaxFunEvals', 2000000);
    %'TolCon', .0001, ...
    %'TolX', .0001);

itersTot=0;
exitflag = 0;
while exitflag == 0     

    [x_opt, fval, exitflag, output, lambda] = fmincon('J', x, ...
        [], [], [], [], vlb, vub, 'NONLCON', options);
    
    itersTot = itersTot+output.iterations;
    disp(sprintf('Saving after %d iterations', itersTot));
    savname=sprintf('sqp%d.mat', itersTot);
    seconds = toc;
    save(savname, 'x_opt', 'fval', 'exitflag', 'output', 'seconds')
    
    x=x_opt;
    toc
end