function prodObj = mtimes(leftObj, rightObj, varargin)
%MTIMES - Define multiplication of Scalars
%
%   Syntax:
%   C = MTIMES(obj, B) is the Scalar convolution product when B is a Scalar truncated to size of obj.
%
%   C = MTIMES(T, obj) is the Scalar image of T(obj) when T is an Operator (linear operator over sequence space) truncated to size of obj.
%
%   C = MTIMES(a, obj) is the Scalar given by a*obj when a is a double truncated to size of obj.
%
%   C = MTIMES(obj, a) is the Scalar given by a*obj = obj*a when a is a double or intval (intvals require right multiplication) truncated to size of obj.
%
%   C = MTIMES(A,B, 'Full') returns a Scalar with no truncation.
%
%   C = MTIMES(obj, B, 'Recursion') returns on the (m-1)^st coefficient with respect to the first dimension when obj.Truncation(1) = m.
%
%   C = MTIMES(A,B, truncationArray) returns a Scalar truncated in each dimension by truncation Array = [N1,...,Nd]).
%
%   Inputs:
%       leftObj - Scalar, double, or Operator
%       rightObj - Scalar, double, or intval
%       truncationSize - string (Full or Recursion) or integer vector of length d.
%
%   Outputs:
%       productObj - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Mar-2017; Last revision: 23-Aug-2018

if isa(leftObj,'Operator')
    %% left action of Operator on a Scalar
    if ~strcmp(rightObj.Basis, 'Taylor')
        warning('mtimes - only implemented for Taylor basis')
    end
    
    switch size(obj,1)
        case 1
            Lv = leftObj.Matrix*rightObj.col;
            convCoefficient = reshape(Lv,rightObj.Truncation);
        otherwise
            error('Operator matrix and vector of Scalars must have the same size')
    end
    
elseif isa(leftObj,'Scalar') && isa(rightObj,'Scalar')
    %% Multiplication of 2 Scalars
    
    % check inputs and raise exceptions
    if ~strcmp(leftObj.Basis, rightObj.Basis)
        error('mtimes - Scalar multiplication requires same basis')
    end
    
    if ~isequal(leftObj.Dimension, rightObj.Dimension)
        error('mtimes - Scalar multiplication requires Scalars of the same dimension')
    end
    
    convBasis = leftObj.Basis;
    if any(ismember(convBasis, 'Chebyshev'))

        % set truncation mode
        if nargin == 2 && isequal(leftObj.Truncation, rightObj.Truncation) % default behavior is embedding into the same truncation space
            truncMode = 'Embed';
            convTruncation = leftObj.Truncation;
        else
            truncMode = varargin{:}; % full 
            convTruncation = leftObj.Truncation + rightObj.Truncation - 1;
        end
        
        chebyshevIdx = 2; % flip along chebyshev dimenions
        convCoefficient = chebtaylortimes(leftObj.Coefficient, rightObj.Coefficient, chebyshevIdx, truncMode);
        prodObj = Scalar(convCoefficient, convBasis, convTruncation);

    elseif strcmp(convBasis, 'Taylor')
        % set truncation mode
        if nargin == 2 && isequal(leftObj.Truncation, rightObj.Truncation) % default behavior is embedding into the same truncation space
            truncMode = arrayfun(@(dim)1:dim, leftObj.Truncation, 'UniformOutput',false); %{1:N1, 1:N2,..., 1:Nd}
            convTruncation = leftObj.Truncation;
        elseif nargin > 2 && isa(varargin{1}, 'double') % specify truncation indices explicitly as {idx1, idx2,...,idxd}
            [truncMode{1:nargin-2}] = varargin{:};
            convTruncation = cellfun(@(dim)length(dim), varargin);
        else % truncMode specified by string varargin
            truncMode = {varargin{1}}; % when objects have different truncations return full convolution
            convTruncation = leftObj.Truncation + rightObj.Truncation - 1;
        end
        
        % compute coefficients of linear discrete convolution
        if strcmp(leftObj.NumericalClass, 'intval') || strcmp(leftObj.NumericalClass, 'intval')
            % One or both are interval Scalars
            convCoefficient = intvaltimes(leftObj.Coefficient, rightObj.Coefficient, truncMode{:});
        else
            % Both Scalars are double
            convCoefficient = doubletimes(leftObj.Coefficient, rightObj.Coefficient, truncMode{:});
        end
        prodObj = Scalar(convCoefficient, convBasis, convTruncation);
        
    else
        error('This basis is not implemented yet')
    end % if
elseif isa(leftObj, 'Scalar')
    %% Scalar multiplication with field scalar (right)
    prodObj = leftObj;
    prodObj.Coefficient = leftObj.Coefficient*rightObj;
    
elseif isa(rightObj,'Scalar')
    %% Scalar multiplication with field scalar (left)
    prodObj = rightObj;
    prodObj.Coefficient = rightObj.Coefficient*leftObj;
    
end % if
end % mtimes

function chebTaylor = chebtaylortimes(leftCoefficient, rightCoefficient, chebyshevIdx, truncMode)
%CHEBTAYLOR - A hacky bootstrap for using existing convolution for Chebyshev-Taylor convlution. Fix this shit.

% return a copy of the coefficients mirrored along each Chebyshev dimension
mirrorLeftCoefficient = Scalar.flipcoef(leftCoefficient, chebyshevIdx);
mirrorRightCoefficient = Scalar.flipcoef(rightCoefficient, chebyshevIdx);

% THIS ONLY WORKS FOR A SINGLE CHEBYSHEV DIRECTION
leftPaddedCoefficient = cat(2, mirrorLeftCoefficient, leftCoefficient(:,2:end));
rightPaddedCoefficient = cat(2, mirrorRightCoefficient, rightCoefficient(:,2:end));
if isa(leftCoefficient, 'intval') || isa(rightCoefficient, 'intval')
    convFullCoefficient = intvaltimes(leftPaddedCoefficient, rightPaddedCoefficient, 'Full');
else
    convFullCoefficient = doubletimes(leftPaddedCoefficient, rightPaddedCoefficient, 'Full');
end

N = size(leftCoefficient,2);
M = size(rightCoefficient,2);
% first meaningful coefficient in the product is the (M+N-1)^st (thanks Max)
if strcmp(truncMode, 'Full')
    chebTaylor = convFullCoefficient(:, N+M-1:end);
elseif strcmp(truncMode, 'Embed')
    chebTaylor = convFullCoefficient(:, N+M-1:N+M-1+min(M,N));
else
    error('This is not implemented yet')
end
end





% Revision History:
%{
24-Jun-2017 - support for 1-dimensional Scalars
02-Jul-2017 - change multiplication to FFT based convolution
19-Aug-2017 - reverted back to non-FFT multiplication. Rigorous and numeric FFT needs to be impleneted more carefully before it can be trusted.
13-Aug-2018 - updated for Scalar class
23-Aug-2018 - Complete overhaul of this method including:
    separation of routine for interval and double. Both routines now return arrays instead of scalars. Assembly of Scalar products is performed
        inside this function. This allows the fast conv routines to be used recursively without calling Scalar methods (much faster).
    implementation for arbitrary dimensions
    support for arbitrary specified truncation
    efficient circular conv routine with padding makes returning smaller truncations much faster
22-May-2019 - Added support for Chebyshev multiplication. It needs work.
%}
