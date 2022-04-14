function evalObj = eval(obj, data)
%EVAL - evaluates the Scalar on specified data.
%
%   EVAL(obj, data) = f(data) where f is the function defined by the Scalar and data is an array of samples on D^d.
%
%   Inputs:
%       obj - Scalar object of dimension d.
%       data - nSamples-by-d array: double or intval evaluation nodes.
%
%   Outputs:
%       evalObj - nSamples length vector: double or intval
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 18-Aug-2018

if ~strcmp(obj(1).Basis, 'Taylor')
    warning('eval - only implemented for Taylor basis')
end


if length(obj) > 1 % vectorized method
    evalObj = cell(size(obj));
    for j = 1:length(obj)
        evalObj{j} = obj(j).eval(data);
    end

else
    switch obj.Dimension
        case 0 % constant function
            evalObj = obj.Coefficient(1)*ones(size(data));
        case 1 % data reshaped to m length row vector
            s = reshape(data,1,[]); % evaluation points as a row vector
            coefficient = reshape(obj.Coefficient,1,[]); % coefficients as a row vector
            S = bsxfun(@power,s,(0:obj.Truncation-1)'); % The jth column of S are the increasing powers of s_j
            evalObj = coefficient*S; % left multiply by row of coefficients
            
        case 2 % data should be an m-by-2 array
            s = data(:,1);
            t = data(:,2);
            S = bsxfun(@power,s',(0:obj.Truncation(2)-1)');
            T = bsxfun(@power,t,0:obj.Truncation(1)-1);
            PS = obj.Coefficient*S;
            TPS = T'.*PS;
            evalObj = sum(TPS,1)';
        otherwise
            error('eval not implemented for this dimension')
    end
end
end %  eval

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
18-Aug-2018 - changed 1-d formula to be consistent with column-major Scalar change
%}
