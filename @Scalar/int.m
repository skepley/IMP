function intObj = int(obj, varargin)
%INT - shift Scalar coefficients by one in the first dimension.
%
%   Syntax:
%       intObj = INT(obj) returns a Scalar which has been "integrated" in the first dimension. On the coefficient level, the m^th coefficient
%       has been shifted and rescaled by 1/m. This is embedded into the truncation space which has 1 additional coefficient and the first
%       coefficient is padded by zero.
%
%       intObj = INT(obj, n) integrates the Scalar n times. 
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       intObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if ~strcmp(obj.Basis, 'Taylor')
    error('int - non-Taylor bases not implemented')
end

if nargin > 1
    nIntegration = varargin{1}; % iterate operator n times
else
    nIntegration = 1;
end

if nIntegration > 1
    obj = int(obj, nIntegration-1);
end

M = obj.Truncation(1);
if strcmp(obj.NumericalClass, 'intval')
    scaleBy = intval(repmat((1:M)', [1, obj.Truncation(2:end)]));
    scaleCoefficient = obj.Coefficient./scaleBy;
    shiftCoefficient = cat(1, intval(zeros([1, obj.Truncation(2:end)])), scaleCoefficient);
elseif strcmp(obj.NumericalClass, 'double')
    scaleBy = repmat((1:M)', [1, obj.Truncation(2:end)]);
    scaleCoefficient = obj.Coefficient./scaleBy;
    shiftCoefficient = cat(1, zeros([1, obj.Truncation(2:end)]), scaleCoefficient);
end
intObj = Scalar(shiftCoefficient, obj.Basis);
end %  shift

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
18-Aug-2018 - changed name from "shift" to "int", added coefficient rescaling to the old shifting operation
20-Aug-2018 - support for arbitrary dimensions
%}
