function reObj = real(obj)
%REAL - return the real part of a Scalar.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

reObj = Scalar(real(obj.Coefficient), obj.Basis);
end %  real

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
