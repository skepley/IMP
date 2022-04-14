function coefDecay = decay(obj)
%DECAY - returns the vector of coefficient norms in the first dimension
%
%   coefDecay = DECAY(oby) returns the ell_1 norm of each coefficient of the Scalar with respect to the first dimension.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       coefDecay - double
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 09-Jun-2017; Last revision: 16-Jan-2019

%   TODO:
%   Add support for Fourier and Chebyshev coefficients.

if numel(obj) > 1 % vectorized method
    if isrow(obj)
       error('decay output only tested for columns of Scalars') 
    end
    rowDecay = arrayfun(@(j)obj(j).decay, 1:numel(obj), 'UniformOutput', false);
    coefDecay = cell2mat(rowDecay);
    
elseif strcmp(obj.NumericalClass, 'Scalar')
    coefDecay = obj.Coefficient(end).norm();
    
else % single Scalar
    coefDecay = sumdim(mid(obj.Coefficient), 2:obj.Dimension);
end
end % decay

function sumDim = sumdim(coef, dimension)
% sum across each dimension specified

if isequal(length(dimension), 0) % return absolute value of vector
    sumDim = abs(coef);
elseif isequal(length(dimension), 1) % sum specified dimension and return
    sumDim = sum(abs(coef), dimension);
else
    sumDim = sumdim(sum(abs(coef),dimension(end)), dimension(1:end-1)); % sum last dimension and recurse
end
end

% Revision History:
%{
02-Jul-2017 - support for Scalar coefficient type
13-Aug-2018 - updated for Scalar class
21-Aug-2018 - support for arbitrary dimensions
16-Jan-2019 - fixed the vectorization. Required removing support for vectorization on row vectors of Scalars.
23-May-2019 - Added support for weighted ell_1 norms so this function works with Fourier or Chebyshev coefficients. 
%}

