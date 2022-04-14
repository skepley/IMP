function convCoefficient = intvaltimes(leftFactor, rightFactor, varargin)
%INTVALTIMES - Scalar multiplication (linear convolution) for intval coefficients
%
%   Description:
%       convObj = INTVALTIMES(leftFactor, rightObj) - description
%
%   Inputs:
%       leftFactor - Scalar leftFactorect with intval coefficients
%       rightObj - Scalar leftFactorect with intval coefficients
%
%   Outputs:
%       convObj - Scalar corresponding to linear convolution leftFactor*rightObj
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 25-Jul-2018; Last revision: 23-Aug-2018

%   TODO:
%   support for other bases

% unpack truncation settings
if nargin == 2 || strcmp(varargin{1}, 'Full')
    truncMode = 'Full';
else
    truncMode = 'Partial';
    [truncSize{1:nargin-2}] = varargin{:};
end

% convert polynomial specified as a coefficient array to coefficient/exponent vectors
dimension = Scalar.dims(leftFactor);
[leftCoefficient, leftExponent] = coef2exponent(leftFactor, dimension);
[rightCoefficient, rightExponent] = coef2exponent(rightFactor, dimension);

% compute discrete linear convolution
switch dimension
    case 0 % multiplication of constants
        convCoefficient = leftCoefficient*rightCoefficient;
        
    case 1 % 1-d discrete convolution
        isColumn = size(leftFactor,2) == 1;
        leftExpArray = repmat(leftExponent, [1 length(rightCoefficient)]);
        rightExpArray = repmat(rightExponent', [length(leftCoefficient),1]);
        
        % build product exponent and coefficients
        convExponent = 1 + leftExpArray + rightExpArray;
        coefArray = leftCoefficient*rightCoefficient.'; % outer product
        convCoefficient = full(sparse(convExponent, 1, coefArray));
        
        % reshape to input size
        if isColumn
            convCoefficient = reshape(convCoefficient, [], 1);
        else
            convCoefficient = reshape(convCoefficient, 1, []);
        end
        
    otherwise % n-d discrete convolution
        leftExpArray = repmat(leftExponent, [1, 1, length(rightCoefficient)]); %1-by-1-by-#Coefficients in right factor
        rightExpArray = shiftdim(repmat(rightExponent', [1, 1, length(leftCoefficient)]), 2); %1-by-1-by-#Coefficients in left factor
        
        % build product exponent and coefficients
        convExponent = reshape(shiftdim(leftExpArray + rightExpArray, 2), [], dimension); % vectorized outer sum of exponent vectors
        outerProduct = reshape(rightCoefficient*leftCoefficient.', [], 1); % vectorized outer product of coefficient vectors 

        base = 1 + max(convExponent, [], 1);
        radix = [1, base(1:end-1)];
        radixIndex = cumprod(radix);
        linearIdx = 1 + convExponent*radixIndex';
        coefVector = sparse(linearIdx, 1 , outerProduct) ; % collect like terms
        convCoefficient = reshape(full(coefVector), base);        
end

% truncate resulting convolution
if strcmp(truncMode, 'Partial')
    convCoefficient = convCoefficient(truncSize{:});
end
end % intvaltimes


function [coef, expon] = coef2exponent(coefArray, dimension)
% converts a coefficient array to a structure of coefficients and exponents

if isequal(dimension, 0)
    coef = coefArray;
    expon = 0;
else
    dims = size(coefArray);
    truncation = dims(dims > 1);
    coef = reshape(coefArray, [], 1); % coefficient as column vector
    % exponent array
    exponVector = arrayfun(@(dim)0:dim, (truncation - 1), 'UniformOutput', false);
    [exponCell{1:dimension}] = ndgrid(exponVector{:});
    for iDim = 1:dimension
        expon(:, iDim) = reshape(exponCell{iDim}, [], 1);
    end
end
end

% Revision History:
%{
13-Aug-2018 - updated for Scalar class
18-Aug-2018 - support for 1-d convolution
22-Aug-2018 - support for arbitrary dimensions
    removed from mtimes method. Functionality now takes intval array as input and output. This function is called inside SCalar/mtimes 
    but may be called outside as well.
19-Mar-2019 - Fixed a bug where complex Cauchy products were wrong. 
%}
