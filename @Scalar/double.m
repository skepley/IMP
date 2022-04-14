function doubleObj = double(obj)
%DOUBLE - convert a Scalar to have double coefficients.
%
%   Syntax:
%       doubleObj = DOUBLE(obj) returns a Scalar with double coefficients from a scalar with any coefficients.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       doubleObj - Scalar with double coefficients
%
%   Subfunctions: none
%   Classes required: Scalar
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018



switch obj.NumericalClass

case 'intval'
  if length(obj) > 1 % vectorized method
      doubleObj = repmat(Scalar(0, obj.Basis, obj(1).Truncation), size(obj));
      for j = 1:length(obj)
          doubleObj(j) = double(obj(j));
      end

  else % convert a single Scalar
      doubleObj = Scalar(mid(obj.Coefficient), obj.Basis);
  end

case 'double'
  doubleObj = obj;

case 'Scalar'
  error('Scalar to double conversion is not yet implemented')

end %  switch
end %  double



% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
