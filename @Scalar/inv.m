function invObj = inv(obj)
%INV - Compute inverse of a Scalar by automatic differentiation
%
%   invObj = INV(obj) returns a Scalar, B, which satisfies B*obj = 1. This is carried out by automatic differentiation by solving
%   the differential equation u' = -obj'*u^2. invObj has same size and numericalClass.
%
%   Inputs:
%       obj - Scalar
%
%   Outputs:
%       invObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-May-2018; Last revision: 27-Mar-2019

%   TODO:
%   Add validation and make the argument recursive.
%   Implement for Fourier/Chebyshev bases.

if ~strcmp(obj.Basis, 'Taylor')
    warning('inv - only implemented for Taylor basis')
end


if strcmp(obj.NumericalClass,'Scalar')
    const = obj.Coefficient(1).Coefficient(1);
else
    const = obj.Coefficient(1);
end

if abs(const) < 1e-13 % check if constant term is near 0.
    disp('Scalar is not invertible. Inverse should not be trusted')
end

switch obj.Dimension
    case 0 % inverse of a real scalar
        a0 = obj.Coefficient(1);
        invObj = Scalar(1/a0, obj.Basis);
        
    case 1 % inverse of 1d power series
        ddtObj = obj.dt;
        c = ddtObj.Coefficient; % expansion for a'
        a0 = obj.Coefficient(1);
        invObj = Scalar(1/a0, obj.Basis, obj.Truncation);
        invObj.Coefficient = reshape(invObj.Coefficient, size(obj.Coefficient)); % match row or column layout of input
        for m = 1:obj.Truncation - 1
            uu = invObj*invObj;
            aPrime = c(m:-1:1);
            invNext = -(1/m)*dot(uu.Coefficient(1:m),aPrime);
            invObj.Coefficient(m+1) = invNext;
        end
        
    case 2 % inverse of 2d power series
        % deg = size(obj.Coefficient);
        ddtObj = obj.dt;
        N = obj.Truncation(2);
        % ddtObj.Coefficient = [ddtObj.Coefficient(2:end,:);zeros(1,deg(2))]; % This padding shouldn't be necessary.
        a0 = Scalar(obj.Coefficient(1,:), obj.Basis);
        invObj = Scalar(inv(a0), obj.Basis, [1, a0.Truncation]);
        for m = 1:obj.Truncation(1)-1
            uu = invObj*invObj;
            aPrime = Scalar(ddtObj.Coefficient(1:m,:), obj.Basis, uu.Truncation);
            % uu_aPrime = uu*aPrime;
            % newCoefficient = -(1/m)*mtimes(uu,aPrime,'Recursion');
            % newCoefficient = -(1/m)*uu_aPrime.Coefficient(m,:);
            invNext =  -(1/m)*mtimes(uu,aPrime,m,1:N);
            invObj.Coefficient(end+1,:) = invNext.Coefficient; % append next coefficient
            invObj.Truncation(1) = invObj.Truncation(1) + 1; % increment truncation dimension
            % invObj.append(newCoefficient);
        end
    otherwise
        error('not yet implemented')
end
end % inv

% Revision History:
%{
02-Jun-2017 - support for 2-dimensional Scalars
09-Jul-2017 - support for Scalar coefficient type
13-Aug-2018 - updated for Scalar class
27-Mar-2019 - Fixed implementation for new Scalar class. 
%}

