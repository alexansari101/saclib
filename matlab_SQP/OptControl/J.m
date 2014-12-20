function [J,g] = J(u)
% Objective function for minimization

% constants
dt = 0.01;
tf = 3.6;
Q = [1000, 0; 0, 10];
R = [0.3];
ulen = length(u);

x=intfn(u);

% takes a state vector, x = [ theta(t_i), theta'(t_i), u(t_i) ];
l=@(x)(([sin(x(1)/2.0), x(2)]-[0, 0])*Q*([sin(x(1)/2.0); x(2)]-[0; 0]) + x(3)*R*x(3));

J = sum(cell2mat( arrayfun(@(x1,x2,u1) l([x1,x2,u1]), x(1:ulen,1),...
    x(1:ulen,2), u(1:ulen)', 'UniformOutput', false) ))/2.0;
% J = cell2mat( arrayfun(@(x1,x2,u1) l([x1,x2,u1]), x(1:ulen,1),...
%     x(1:ulen,2), u(1:ulen), 'UniformOutput', false) );

if nargout > 1 % fun called with two output arguments
 g =  3*u'/10; % Compute the gradient evaluated at x
end
