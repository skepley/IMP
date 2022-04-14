function varargout = sqrt(obj, varargin)
%SQRT - compute square roots for Scalars using automatic differentiation.
%
%   B = SQRT(obj) is a scalar with the same truncation embedding as obj which approximately satisfies B*B = obj. The constant term for B is positive.
%
%   B = SQRT(obj, -1) is the negative solution branch (i.e. B*B = obj and the constant term for B is negative).
%
%   B = SQRT(obj) is a Scalar with the same truncation embedding as obj which approximately satisfies B*B = obj.
%
%   [B, inv(B)] is a pair of Scalars satisfying B*B = obj and B*inv(B) = 1.
%
%   Inputs:
%       obj - Scalar
%       branch - integer either {1, -1}
%
%   Outputs:
%       sqrtObj - Scalar
%       invObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-May-2018; Last revision: 13-Aug-2018

if ~strcmp(obj.Basis, 'Taylor')
    error('sqrt - only implemented for Taylor basis')
end

if strcmp(obj.NumericalClass,'intval')
    warning('sqrt - untested for interval coefficients')
end

if abs(obj.Coefficient(1,1)) < 1e-9
    warning('constant term < 1e-9. square root should not be trusted')
end

if nargin ==2
    branch = varargin{1};
else
    branch = 1;
end

switch obj.Dimension
    case 0 % sqrt of a field scalar
        if isreal(obj.Coefficient(1)) % choose specified branch of real square root
            x0 = abs(obj.Coefficient(1));
            sqrtObj = Scalar(branch*sqrt(x0), obj.Basis);
            invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis);
        else
            x0 = obj.Coefficient(1); % CHOOSES THE PRINCIPAL BRANCH OF THE COMPLEX SQUARE ROOT
            sqrtObj = Scalar(sqrt(x0), obj.Basis);
            invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis);
            warning('Complex valued Scalars always return the principal branch')
        end
        
    case 1 % sqrt of 1d power series
        xDot = obj.dt;
        
        if isreal(obj.Coefficient(1))
            x0 = abs(obj.Coefficient(1)); % |x_0|
            sqrtObj = Scalar(branch*sqrt(x0), obj.Basis, 1); % y
            invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis, 1); % z
            % sqrtObj = Scalar(branch*sqrt(x0), obj.Basis, obj.Truncation); % y
            % invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis, obj.Truncation); % z
            % invObj.Coefficient = reshape(invObj.Coefficient, [], 1); % make into column layout
        else
            x0 = obj.Coefficient(1);
            sqrtObj = Scalar(sqrt(x0), obj.Basis, 1);
            invObj = Scalar(1/sqrtObj.Coefficient(1), obj.Basis, 1);
        end
        
        % recursively solve differential equation for the coefficients
        for m = 1:obj.Truncation - 1
            xDotTruncation = Scalar(xDot.Coefficient(1:m), obj.Basis, invObj.Truncation);
            % xDotTruncation = Scalar(xDot.Coefficient(1:m), obj.Basis);
            
            dsqrt = xDotTruncation*invObj; % xDot*z
            zz = invObj*invObj; % z*z
            sqrtNext = (0.5/m)*dsqrt.Coefficient(m); % this is a double
            invNext = (-0.5/m)*mtimes(zz, dsqrt,m); % this is a Scalar
            
            % append next coefficient and increment truncation dimension
            invObj.Coefficient(end+1) = invNext.Coefficient;
            sqrtObj.Coefficient(end+1) = sqrtNext;
            invObj.Truncation(1) = invObj.Truncation(1) + 1;
            sqrtObj.Truncation(1) = sqrtObj.Truncation(1) + 1;
        end
        % reshape to row or column layout to match input
        invObj.Coefficient = reshape(invObj.Coefficient, size(obj.Coefficient));
        sqrtObj.Coefficient = reshape(sqrtObj.Coefficient, size(obj.Coefficient));
        
    case 2 % sqrt of 2d power series
        xDot = obj.dt;
        N = obj.Truncation(2);
        x0 = Scalar(obj.Coefficient(1,:), obj.Basis);
        [y0,z0] = sqrt(x0);
        sqrtObj = Scalar(y0, obj.Basis, [1, x0.Truncation]);
        invObj = Scalar(z0, obj.Basis, [1, x0.Truncation]);
        for m = 1:obj.Truncation(1)-1
            xDotTruncation = Scalar(xDot.Coefficient(1:m,:), obj.Basis, invObj.Truncation);
            dsqrt = xDotTruncation*invObj; % xDot*z
            zz = invObj*invObj; % z*z
            sqrtNext = (0.5/m)*dsqrt.Coefficient(m,:); % this is a double
            invNext = (-0.5/m)*mtimes(zz, dsqrt, m, 1:N); % this is a Scalar
            
            % append next coefficient and increment truncation dimension
            invObj.Coefficient(end+1,:) = invNext.Coefficient;
            sqrtObj.Coefficient(end+1,:) = sqrtNext;
            invObj.Truncation(1) = invObj.Truncation(1) + 1;
            sqrtObj.Truncation(1) = sqrtObj.Truncation(1) + 1;
        end
    otherwise
        error('sqrt not implemented for this dimension')
end

if nargout == 1
    varargout{1} = sqrtObj;
elseif nargout ==2
    varargout{1} = sqrtObj;
    varargout{2} = invObj;
end
end

% Revision History:
%{
17-Jun-2017 - computation rewritten to use automatic differentiation
    support for lower dimensional square roots
    added option to return inv(obj) which is computed as a by product of automatic differentiation
13-Aug-2018 - updated for Scalar class
Marchish? 2019 - Fixed several bugs with the new Scalar class.
19-May-2019 - Added support for complex valued Scalar square roots. This only chooses the principal branch.
%}
