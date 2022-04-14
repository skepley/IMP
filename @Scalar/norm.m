function scalarNorm = norm(obj, varargin)
%NORM - computes the weighted ell_1 norm of a Scalar
%
%   NORM returns the weighted ell_1 norm for a Scalar along specified dimensions. The weights are specified in the Scalar properties.
%
%   Syntax:
%       scalarNorm = NORM(obj) returns the ell_1 norm of obj
%       scalarNorm = NORM(obj, dims) returns the ell_1 norm of obj in the directions specified by dims
%
%   Inputs:
%       obj - Scalar
%       dims - Vector of dimensions along which to sum
%
%   Outputs:
%       scalarNorm - double or intval depending on Scalar NumericalClass
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

% parse input
p = inputParser;
addOptional(p, 'Weight', ones(1, obj.Dimension))
addParameter(p,'Dimension', 'All')

parse(p, varargin{:})
ellOneWeight = p.Results.Weight; % A vector of ell_1 weights: [nu_1,...,nu_d] or 'ones' when nu_i = 1 for all i
dimension = p.Results.Dimension;

if numel(obj) > 1 % vectorized norm
    for j = 1:numel(obj)
        scalarNorm(j) = obj(j).norm(varargin{:});
    end
    scalarNorm = reshape(scalarNorm,size(obj));
    
else % norm of a single scalar
    if all(ellOneWeight == 1) % All weights equal to 1
        if strcmp(dimension, 'All')
            scalarNorm = sum(abs(obj.Coefficient(:)));
        else
            scalarNorm = sum(abs(obj.Coefficient), dimension);
        end
        
    elseif strcmp(dimension, 'All')
        partialNorm = abs(obj.Coefficient);
        for wIdx = 1:obj.Dimension
            weight = ellOneWeight(wIdx);
            partialNorm = collapseDimension(weight, partialNorm, strcmp(obj.Basis{wIdx}, 'Chebyshev'));
        end
        scalarNorm = partialNorm;
    else
        error('This is not implemented yet')
    end % if
end % if
end % norm

function partialNorm = collapseDimension(weight, partialNorm, isCheb)
% collapse by evaluation in the first dimension

sz = size(partialNorm);
if Scalar.dims(partialNorm) == 1
    sz = sz(sz > 1); % remove the singleton dimension
    partialNorm = partialNorm.';
end
N = sz(1);

if isequal(weight, 1)
    partialNorm = sum(partialNorm, 1);
    
    % get a matrix of weights to multiply in the first dimension
elseif isCheb % multiply all coefficients after the constant by a factor of 2
    W = cat(1, ones([1, sz(2:end)]), 2*repmat(weight.^(1:N-1).', [1, sz(2:end)]));
    partialNorm = squeeze(sum(W.*partialNorm, 1)); % evaluate in the dimension 1
    
else
    W = repmat(weight.^(0:N-1).', [1, sz(2:end)]);
    partialNorm = squeeze(sum(W.*partialNorm, 1)); % evaluate in the dimension 1
    
end

end
% Revision History:
%{
08-Aug-2018 - moved out of classdef file.
22-May-2019 - Added support for arbitrary weights and Chebyshev basis Scalars. Changed this to take weights as an input instead of as a property of Scalar.
%}
