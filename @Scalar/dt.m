function dObj_dt = dt(obj)
%DT - compute the derivative of a Scalar
%
%   DT(obj) returns a Scalar corresponding to the derivative of obj in the first dimension. The derivative is embedded in the same
%   truncation space by padding with zeros in the last coordinate.
%
%   Inputs:
%       obj - Scalar object
%
%   Outputs:
%       dObj_dt - Scalar object
%
%   Subfunctions: none
%   Classes required: Scalar
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 29-Mar-2019
%   
% TODO: 
%       1. Update this function for the new Scalar class including support for any dimension
%       2. Implement support for taking derivatives along arbitrary dimensions instead of just the first. 

% warning('this version of dt includes the shift. check the code which called this')
if obj.Dimension ==1 % space is 1-dim
    C = reshape(0:obj.Truncation(1)-1, size(obj.Coefficient));
    scaleCoefficient = C.*obj.Coefficient; % includes the leading zero
    scaleCoefficient(end+1) = 0; % pad last coefficient
    dObj_dt = Scalar(scaleCoefficient(2:end), obj.Basis); % shift and return Scalar

elseif obj.Dimension == 2 % space is 1-dim
    C = repmat((0:obj.Truncation(1)-1),obj.Truncation(2),1)';
    scaleCoefficient = C.*obj.Coefficient; % includes the leading zero
    scaleCoefficient(end+1,:) = 0; % pad last coefficient with a row of zeros
    dObj_dt = Scalar(scaleCoefficient(2:end,:), obj.Basis); % shift and return Scalar

elseif obj.Dimension ==3 % space is 2-dim
    C = zeros(obj.Truncation);
    for j = 2:obj.Truncation(1)
        C(j,:,:) = (j-1)*ones(obj.Truncation(2:3));
    end
    dObj_dt = Scalar(C.*obj.Coefficient, obj.Basis);
end
end %  dt

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
29-Mar-2019 - Fixed a bug for 1-dimensional Scalars which assumed Scalar coefficients were in row layout. Now the coefficients are always reshaped to match either. 
    Fixed another bug where coefficients were being scaled but not shifted. 
%}
