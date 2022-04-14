function operator = leftmultiplicationoperator(obj)
%LEFTMULTIPLICATIONOPERATOR - returns multiplication operator associated with a Scalar.
%
%   LEFTMULTIPLICATIONOPERATOR(obj) returns a linear operator, A, whose action on any Scalar, h, is given by convolution with obj.
%   In other words, for every h, A(h) = obj*h.
%
%   Inputs:
%       obj - Scalar object
%
%   Outputs:
%       operator - Operator object
%
%   Subfunctions: none
%   Classes required: Operator
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if ~strcmp(obj.Basis, 'Taylor')
    warning('leftmultiplicationoperator - only implemented for Taylor basis')
end


if strcmp(obj.NumericalClass, 'intval')
    warning('this method is untested with interval valued Scalars')
end


padblock = [zeros(obj.Truncation(1)-1, obj.Truncation(2)-1), zeros(obj.Truncation(1)-1,obj.Truncation(2));...
    zeros(obj.Truncation(1), obj.Truncation(2)-1),obj.Coefficient];

leftMultiplier(obj.Truncation(1), obj.Truncation(2)) = Scalar(0, obj.Basis, obj.Truncation);
for j = 1:obj.Truncation(1)
    for k = 1:obj.Truncation(2)
        leftMultiplier(obj.Truncation(1)+1-j,obj.Truncation(2)+1-k) = Scalar(padblock(j:j+obj.Truncation(1)-1, k:k+obj.Truncation(2)-1), obj.Basis);
    end
end
operator = Operator(leftMultiplier, obj.Truncation);
end % leftmultiplicationoperator

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
