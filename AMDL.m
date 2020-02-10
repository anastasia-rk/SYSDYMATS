function [result] = AMDL(r_n,nData,nTerms)
% Notation mataches that of (Wei,Lang,Billings, 2008) and differs from the main file
% Input
%   nData  - length of the dataset;
%   nTerms - number of model terms;
%   r_n    - associated model residual;
% Output
%   result - AMDL associated with the k-th dataset
%   AMDL - Approximate Minimum Description Length criterion

    norm_r = norm(r_n);
    %     result = 2*log(nData)*nTerms/nData + log(norm_r/nData);
     result  = 1.5*log(nData)*nTerms/nData + 0.5*log(norm_r/nData);
end