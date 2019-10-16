function [result] = cor_sqr(x,y)
% Notation mataches that of (Wei,Lang,Billings, 2008) and differs from the main file
% Input
%   x - vector 1 of size n x 1
%   y - vector 2 of size n x 1
% Output
fact1 = (x'*y)^2;
fact2 = (x'*x);
fact3 = (y'*y);
result = fact1/(fact2*fact3);
end