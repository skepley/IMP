function polyObj = intlabpoly(obj, varargin)
%INTLABPOLY - % returns intlab polynomial object with intval coefficients corresponding to Scalar.
%
%   Syntax:
%       polyObj = INTLABPOLY(obj) returns an IntLab polynom object for the given Scalar obj.
%       polyObj = INTLABPOLY(obj, vars) returns an IntLab polynom with variable names specified as strings in vars.
%
%   Inputs:
%       obj - Scalar
%       vars - cell vector of strings
%
%   Outputs:
%       polyObj - IntLab polynom
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: IntLab toolbox

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018



if numel(obj) > 1 % vectorize function call
    polyObj = arrayfun(@(scalar)intlabpoly(scalar,varargin{:}), obj, 'UniformOutput',false);

else % convert a single Scalar
    if strcmp(obj.NumericalClass, 'Scalar')
        warning('Scalar NumericalClass is deprecated')
        homogObj = homog(obj);
    else
        homogObj = obj;
    end

    switch obj.Dimension
        case 1 % NumericalClass can only be double or intval
            exponObj = (0:obj.Truncation - 1)';
            if strcmp(homogObj.NumericalClass,'double')
                coefObj = reshape(midrad(homogObj.Coefficient,0),[],1);
            elseif strcmp(homogObj.NumericalClass,'intval')
                coefObj = reshape(homogObj.Coefficient,[],1);
            end
            polyObj = polynom(coefObj,exponObj,'s');
        case 2
            if nargin > 1
                vars = varargin{1};
            else
                vars = {'s','t'};
            end

            [S,T] = meshgrid(0:homogObj.Truncation(2)-1,0:homogObj.Truncation(1)-1);
            exponObj = [reshape(S,[],1),reshape(T,[],1)];

            if strcmp(homogObj.NumericalClass,'double')
                coefObj = reshape(midrad(homogObj.Coefficient,0),[],1);
            elseif strcmp(homogObj.NumericalClass,'intval')
                coefObj = reshape(homogObj.Coefficient,[],1);
            end
            polyObj = polynom(coefObj,exponObj,vars);
        otherwise
          error('intlabpoly not yet implemented for Dimension other than 1')
    end % switch
end % if
end % intlabpoly

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
