function intvalObj = intval(obj)
%INTVAL - convert a Scalar to have intval coefficients.
%
%   Syntax:
%       intvalObj = INTVAL(obj) returns a Scalar with intval coefficients from a scalar with any coefficients.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       intvalObj - Scalar with intval coefficients
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: IntLab

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

switch obj.NumericalClass
  case 'double'
    if length(obj) > 1 % vectorized method
        intvalObj = repmat(Scalar(0, obj.Basis, obj(1).Truncation), size(obj));
        for j = 1:length(obj)
            intvalObj(j) = obj(j).intval;
        end
    else % convert a single Scalar
        intvalObj = Scalar(midrad(obj.Coefficient, 0), obj.Basis);

        % replaced midrad(0,0) with midrad(0,eps) so that intLab will keep placeholder zeros
        intvalObj.Coefficient(obj.Coefficient == 0) = midrad(0, eps);
    end
  case 'intval'
    intvalObj = obj;
  case 'Scalar'
    error('Scalar to intval conversion is not yet implemented')
end %  switch
end %  intval

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
