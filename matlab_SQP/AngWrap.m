function [x] = AngWrap(x)
% x is a vector [theta, theta'...] whose first element is an angle
% xwrap is a wrapped vector where theta in [-pi, pi)
% angles are asssumed to be in radians

x(1) = mod(x(1)+pi, 2*pi);

if x(1) < 0
    x(1) = x(1) + 2*pi;
end

x(1) = x(1) - pi;