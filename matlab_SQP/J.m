function [J,g] = J(x)
% Objective function for minimization

global dt xlen

% constants
Q = [1000, 0; 0, 10];
R = [0.3];
xdes = [0, 0];

%%% NON - ANGLE WRAPPED VERSION 1 %%%
% 
% % takes a state vector, x = [ theta(t_i), theta'(t_i), u(t_i) ];
% l=@(x)(([sin(x(1)/2.0), x(2)]-xdes)*Q*([sin(x(1)/2.0); x(2)]-xdes') ... 
%     + x(3)*R*x(3));
% 
% J = dt*sum(cell2mat( arrayfun(@(x1,x2,u1) l([x1,x2,u1]), x(1:3:xlen),...
%     x(2:3:xlen), x(3:3:xlen), 'UniformOutput', false) ))/2.0;
% 
% % Compute the gradient evaluated at x
% if nargout > 1 % fun called with two output arguments
%  g(1:3:xlen) = dt*Q(1,1)*sin(x(1:3:end))/4.0;
%  g(2:3:xlen) = dt*Q(2,2)*x(2:3:end);
%  g(3:3:xlen) = dt*R(1,1)*x(3:3:end); 
%  % g = [ dt*Q(1,1)*sin(x(1,:))/4.0; dt*Q(2,2)*x(2,:); dt*R(1,1)*x(3,:) ];
% end

%%% FOR REFERENCE
%%% A = [0, 1; (g*cos(x(1)) - sin(x(1))*u(1))/h, 0]
%%% B = [0; cos(x(1))/h];
%%%

%%% ANGLE WRAPPED VERSION 2 %%%
%
x(1:3:end)=cell2mat( arrayfun(@(x1) AngWrap([x1]), x(1:3:end),...
    'UniformOutput', false) );

% takes a state vector, x = [ theta(t_i), theta'(t_i), u(t_i) ];
l=@(x)(([x(1), x(2)]-xdes)*Q*([x(1); x(2)]-xdes') ... 
    + x(3)*R*x(3));

J = dt*sum(cell2mat( arrayfun(@(x1,x2,u1) l([x1,x2,u1]), x(1:3:xlen),...
    x(2:3:xlen), x(3:3:xlen), 'UniformOutput', false) ))/2.0;

% Compute the gradient evaluated at x
if nargout > 1 % fun called with two output arguments
 g(1:3:xlen) = dt*Q(1,1)*x(1:3:end);
 g(2:3:xlen) = dt*Q(2,2)*x(2:3:end);
 g(3:3:xlen) = dt*R(1,1)*x(3:3:end); 
 % g = [ dt*Q(1,1)*sin(x(1,:))/4.0; dt*Q(2,2)*x(2,:); dt*R(1,1)*x(3,:) ];
end