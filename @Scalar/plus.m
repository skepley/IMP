function sumObj = plus(leftObj, rightObj)
%PLUS - Define addition for Scalars
%
%   Syntax:
%   C = PLUS(L, R) when L and R are one Scalar and one field scalar (intval or double). Returns the sum of a Scalar obj with
%     the Scalar embedding of a (double or intval) as a constant coefficient sequence.
%
%   C = PLUS(L, R) when L and R are both Scalars with the same truncation embedding. Returns the sum of the Scalars with
%     the same truncation.
%
%   Inputs:
%       leftObj - Scalar, double, or intval
%       rightObj - Scalar, double, or intval
%
%   Outputs:
%       sumObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-May-2017; Last revision: 08-Aug-2018

if isa(leftObj,'double') % Sum of Scalar and (double) Fscalar (left or right)
    sumObj = rightObj;
    if isa(sumObj.Coefficient,'Scalar')
        sumObj.Coefficient(1).Coef(1) = sumObj.Coefficient(1).Coef(1) + leftObj;
    else
        sumObj.Coefficient(1) = sumObj.Coefficient(1) + leftObj;
    end
    
elseif isa(rightObj,'double') || isa(rightObj,'intval')  % Sum of Scalar and (double or intval) Fscalar (right only)
    sumObj = leftObj;
    if isa(sumObj.Coefficient,'Scalar')
        sumObj.Coefficient(1).Coef(1) = sumObj.Coefficient(1).Coef(1) + rightObj;
    else
        sumObj.Coefficient(1) = sumObj.Coefficient(1) + rightObj;
    end
    
else % Sum of 2 Scalars
    if ~isequal(leftObj.Dimension, rightObj.Dimension)
        error('Addition for different surface Dimensions is not defined')
    end
    
    if ~strcmp(leftObj.Basis, rightObj.Basis)
        error('mtimes - Scalars must have the same basis type')
    else
        basis = leftObj.Basis;
    end
    
    if ~isa(leftObj.Coefficient,'Scalar') && ~isa(rightObj.Coefficient,'Scalar') % both summands have double or intval Coefs
        if leftObj.Truncation == rightObj.Truncation
            sumObj = Scalar(leftObj.Coefficient + rightObj.Coefficient, basis, leftObj.Truncation);
        else % summands have non-equal size.
            padUp = max([leftObj.Truncation; rightObj.Truncation]);
            sumObj = Scalar(leftObj.Coefficient, basis, padUp) + Scalar(rightObj.Coefficient, basis, padUp);
        end
        
    elseif isa(leftObj.Coefficient,'Scalar') && isa(rightObj.Coefficient,'Scalar') % both summands have Scalar Coefs
        if isequal(leftObj.Truncation,rightObj.Truncation) % summands have same size.
            sumObj = leftObj;
            for j = 1:leftObj.Truncation{1}
                sumObj.Coefficient(j) = sumObj.Coefficient(j) + rightObj.Coefficient(j);
            end
        else % summands have non-equal size.
            commonModes = min([leftObj.Truncation{1};rightObj.Truncation{1}]);
            if leftObj.Truncation{1} > commonModes % left summand is the large one.
                largerSummand = leftObj;
            else
                largerSummand = rightObj;
            end
            
            sumObj = Scalar(largerSummand, basis);
            for j = 1:commonModes
                sumObj.Coefficient(j) = sumObj.Coefficient(j) + rightObj.Coefficient(j); % sum over non-zero modes.
            end
        end
    else
        error('Addition of Scalar CoefType and non-Scalar CoefType not supported')
    end
end
end %  plus

% Revision History:
%{
02-Jul-2017 - support for Scalar coefficients
08-Aug-2018 - moved out of classdef file
%}
