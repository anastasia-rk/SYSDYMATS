function [result] = ERR(x,y)
% Input
%   x - vector 1 of size n x 1
%   y - vector 2 of size n x 1
% Output
fact1 = (x'y)^2;
fact2 = (x'x)^2;
fact3 = (x'y)^2;
result = fact1/(fact2*fact3);
end