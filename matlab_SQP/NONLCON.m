function [C, Ceq, DC, DCeq] = NONLCON(x)
% Returns Nonlinear Inequality and Equality Constraints
% C(x) <= 0  OR  C(x) = 0

global dt xlen h x0 derivC

g = 9.81;

% Inequality Constraints C(x) <= 0
C=[];

% Equality Constraints C(x) = 0
% Ceq=[];

% Compute what the true state should be from initial conditions and
% controls
% u = cat(2,x0(1,3),x(3:3:end)); % u = cat(2,x0(3,1),x(3:3:end));
% xtrue = intfn(x0(1:2),u);
% Ceq = 0.*x;
% Ceq(1:3:end) = x(1:3:end)-xtrue(:,1)';
% Ceq(2:3:end) = x(2:3:end)-xtrue(:,2)'; % constraint satisfied when traj satisfies dynamics
%
% A = [0, 1; (g*cos(x(1)) - sin(x(1))*u(1))/h, 0]
% B = [0; cos(x(1))/h];
%
% Constraint: 0 = x0-x1 + dt*f(x0,u0)
% Ceq(1) = x0(1)-x(1)+dt*(x0(2));
% Ceq(2) = x0(2)-x(2)+dt*( ( g*sin(x0(1)) + x0(3)*cos(x0(1)) ) / h );
% Ceq(3) = 0;
% Ceq(4:3:xlen) = x(1:3:(end-3))-x(4:3:end)+dt*( x(2:3:(end-3)) );
% Ceq(5:3:xlen) = x(2:3:(end-3))-x(5:3:end)+dt*( ( g*sin(x(1:3:(end-3))) ...
%     + x(3:3:(end-3)).*cos(x(1:3:(end-3))) ) / h );
% Ceq(6:3:xlen) = 0;
%
% Only include non-zero constraints
Ceq(1) = x0(1)-x(1)+dt*(x0(2));
Ceq(2) = x0(2)-x(2)+dt*( ( g*sin(x0(1)) + x0(3)*cos(x0(1)) ) / h );
Ceq(3:1+xlen/3) = x(1:3:(end-3))-x(4:3:end)+dt*( x(2:3:(end-3)) );
Ceq(end+1:end+xlen/3-1) = x(2:3:(end-3))-x(5:3:end)+dt*( ( g*sin(x(1:3:(end-3))) ...
    + x(3:3:(end-3)).*cos(x(1:3:(end-3))) ) / h );
%
% Compute the gradient evaluated at x
if nargout > 2 % fun called with 4 output arguments
  % deriv = [ -1 0 0 0 ...; 0 -1 0 0 0 ...; 1 dt 0 -1 0 ...; ... ; 
  %          g*cos(x(1))-sin(x(1))*x(3))/h 1 cos(x(1))/h 0 -1 0 ...;
  %          0 0 0 g*cos(x(4))-sin(x(4))*x(6))/h 1 cos(x(4))/h 0 -1 0 ... ; ]
  % deriv = zeros(length(Ceq), xlen);
  % deriv(1,1) = -1;
  % deriv(2,2) = -1;
  % dc1 = [1 dt 0 -1]; % rows [3 : 1+xlen/3]
  % for i = [3:1+xlen/3]
  %   deriv(i,3*(i-3)+1:3*(i-3)+4) = dc1;
  % end
  % rb = 2+xlen/3;
  % re = length(Ceq);
  % for i = [rb:re]
  %   deriv(i,3*(i-rb)+1:3*(i-rb)+5) = ... 
  %       [(g*cos(x((i-rb)*3+1))-sin(x((i-rb)*3+1))*x((i-rb)*3+3))/h 1 cos(x((i-rb)*3+1))/h 0 -1];
  % end
  % DC = [];
  % DCeq = deriv';
  %
  rb = 2+xlen/3;
  re = length(Ceq);
  for i = [rb:re]
    derivC(i,[3*(i-rb)+1, 3*(i-rb)+3]) = ... 
        [(g*cos(x((i-rb)*3+1))-sin(x((i-rb)*3+1))*x((i-rb)*3+3))/h cos(x((i-rb)*3+1))/h];
  end
  DC = [];
  DCeq = derivC';
end