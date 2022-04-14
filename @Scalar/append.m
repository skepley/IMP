function updateObj = append(obj, nextCoefficient, varargin)
%APPEND - Appends a new coefficient with respect to specified dimension (default is 1st dimension).
%
%   APPEND() - A more detailed description of the function
%
%   Syntax:
%       updateObj = APPEND(obj, nextCoefficient) returns the Scalar object with nextCoefficient appended to the first dimension.
%        nextCoefficient should be obj.Dimension-1 dimensional array of size obj.truncation(2:end).
%
%       updateObj = APPEND(obj, nextCoefficient, dim) returns the Scalar object with nextCoefficient appended to the dimension specified by dim.
%
%   Inputs:
%       obj - Scalar
%       nextCoefficient - Scalar or array of same NumericalClass
%       dim - positive integer
%
%   Outputs:
%       updateObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 09-May-2017; Last revision: 13-Aug-2018

% TODO:
% Fix code for appending along arbitrary dimension

if strcmp(obj.NumericalClass, 'Scalar') % append Scalar coefficient
    obj.Coefficient(end+1) = nextCoefficient;
    obj.Truncation(1) = obj.Truncation(1) + 1;

else % append double/intval array
    if isa(nextCoefficient,'Scalar')
        nextCoefficient = nextCoefficient.Coefficient;
    end

    try
        switch obj.Dimension
            case 1
                obj.Coefficient(end+1, 1) = nextCoefficient;
            case 2
                obj.Coefficient(end+1,:) = nextCoefficient;
            case 3
                obj.Coefficient(end+1,:,:) = nextCoefficient;
            otherwise
                error('Not yet implemented for higher dimension')
        end
        obj.Truncation(1) = obj.Truncation(1) + 1;
    catch
        error('Appended coefficient must have same dimension as the surface being appended to')
    end
end
updateObj = obj;
end %  append

% Revision History:
%{
24-Jun-2017 - support for specifying dimension
02-Jul-2017 - support for Scalar coefficients
13-Aug-2018 - updated for Scalar class
%}
