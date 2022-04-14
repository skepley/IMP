function varargout = exponent(obj)
%EXPONENT - Return a Scalar as an array of coefficients and exponents
%
%   Syntax:
%       output = EXPONENT(input1, input2)
%       output = EXPONENT(input1, input2, input3)
%
%   Description:
%       EXPONENT() - description
%
%   Inputs:
%       obj - instance of Scalar object
%
%   Outputs:
%       coefficient - coefficient array as a column vector
%       exponent - exponent array of size: #length(coefficients)-by-Scalar dimension
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 20-Aug-2018

if isequal(obj.Dimension, 0)
    coefficient = obj.Coefficient;
    exponent = 0;
else
    coefficient = reshape(obj.Coefficient, [], 1); % coefficient as column vector
    % exponent array
    exponVector = arrayfun(@(dim)0:dim-1,obj.Truncation, 'UniformOutput', false);
    exponCell = cell(1, obj.Dimension);
    [exponCell{:}] = ndgrid(exponVector{:});
    for iDim = 1:obj.Dimension
        exponent(:, iDim) = reshape(exponCell{iDim}, [], 1);
    end
end

switch nargout
    case 1
        structOut.Coefficient = coefficient;
        structOut.Exponent = exponent;
        varargout{1} = structOut;
    case 2
        varargout{1} = coefficient;
        varargout{2} = exponent;
end
end %  exponent

% Revision History:
%{
13-Aug-2018 - updated for Scalar class
20-Aug-2018 - support for arbitrary dimensions
%}
