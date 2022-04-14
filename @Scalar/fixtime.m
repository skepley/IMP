function evalObj = fixtime(obj, t)
%FIXTIME - partial evaluation of Scalar along the first (time) dimension.
%
%   evalObj = FIXTIME(obj, tau) returns a Scalar of dimension (d-1) corresponding to evaluation of obj along its first dimension. In usual cases
%       this dimension is utilized by the time variable so that this evaluation is the time-tau map for the Scalar.
%       IMPORTANT NOTE: This evaluation is in MATERIAL time, not real time.
%
%   Inputs:
%       obj - Scalar
%       t - A real number in the interval [-1,1] (for Taylor basis)
%
%   Outputs:
%       evalObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 23-Apr-2017; Last revision: 08-Mar-2018

%   TODO:
%   support for Fourier and Chebyshev bases
%   support for evaluation along other dimensions
%   support for specifying real time evaluation

if ~strcmp(obj(1).Basis, 'Taylor')
    warning('fixtime - only implemented for Taylor basis')
end

if length(obj) ==1 % obj is a single Scalar
    if t==1 % fast evaluation by summing over the first dimension which is equivalent to evaluation at tau = 1.
        evalCoefficient = sum(obj.Coefficient,1);
    elseif t == -1 % fast evaluation by summing over the first dimension with alternating signs which is equivalent to evaluation at tau = -1.
        evalCoefficient = sum(obj.Coefficient(1:2:end),1) - sum(obj.Coefficient(2:2:end),1);
    else
        switch obj.Dimension
            case 1
                timeVector = bsxfun(@power, t, 0:length(obj.Coefficient)-1); % vector of powers
                evalCoefficient = dot(timeVector, obj.Coefficient); % dot product to evaluate
            case 2
                timeVector = bsxfun(@power, t, 0:size(obj.Coefficient,1)-1); % vector of powers along first dimension
                evalCoefficient = timeVector*obj.Coefficient; % left multiply by time_vector to evaluate
            otherwise
                error('Not implemented for arbitrary tau in higher dimensons')
        end
    end
    evalObj = Scalar(squeeze(evalCoefficient), obj.Basis(2:end)); % boundary object is written without the first basis of the interior object
    
else % obj is a vector of Scalars
    evalObj(length(obj)) = obj(length(obj)).fixtime(t); % initialize by fixing time in last coordinate
    for j = 1:length(obj)-1 % loop over first n-1 coordinates fixing time in each
        evalObj(j) = obj(j).fixtime(t);
    end
    evalObj = reshape(evalObj,size(obj)); % change layout to match input
end
end % fixtime

% Revision History:
%{
15-May-2017 - support for 2-dimensional Scalars
13-Aug-2018 - updated for Scalar class
08-Mar-2019 - support for 1-dimensional Scalars and fast evaluation for tau = 1 and -1 cases. Fixed a bug where the output array was not in the same
    layout as the input. 
%}
