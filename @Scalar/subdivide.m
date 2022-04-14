function varargout = subdivide(obj, parmRange, objDomain)
%SUBDIVIDE - rigorously subdivide a Scalar parameterization into a piecewise Scalar parameterization.
%
%   Syntax:
%       subScalar = SUBDIVIDE(obj, parmRange) is a vector of Scalars of length nDivs which is the number of rows of parmRange.
%           Each row of parmRange specifies a subset of D^d on which to parameterize each subScalar.
%
%
%       [subScalar, subDomain] = SUBDIVIDE(obj, parmRange) subDomain is a nDivs-by-(2d) array of material coordinates with respect to [-1,1]^d. These coordinates
%           are consecutively inherited by successive calls of subdivide.
%
%   Inputs:
%       obj - Scalar
%       parmRange - nDivs-by-(2d) array
%       objDomain - A subset of [-1,1]^d on which the Scalar object defines an analytic function
%
%   Outputs:
%       newParm - vector of Scalars
%       newMTC - double array
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jan-2017; Last revision: 8-Mar-2018
%
% TODO: Take advantage of improved error bounds for subScalars in the interior of the spatial domain

if nargin > 3
    error('time-tau evaluation inside Scalar/subdivide was removed. Evaluate this before calling subdivide')
end

scalarDimension = obj(1).Dimension(1); % get dimension of Scalars
scalarBasis = obj(1).Basis; % get Basis for Scalars

if size(parmRange, 2) ~= 2*scalarDimension
    error('parmRange should have 2*dimension-many columns')
end

% use interval subdivision or not
if strcmp(obj(1).NumericalClass,'intval')
    isValid = true;
else
    isValid = false;
end

% Get spatial domain for parent surface
if nargout == 2
    switch scalarDimension
        case 1
            % map parmRange which are computed on [-1,1]^d into the actual domain for these Scalars
            rescale = @(s)(mean(objDomain) + .5*s*diff(objDomain));
            newMTC = rescale(parmRange);
        otherwise
            error('not yet implemented for higher dimensions')
    end
end

nScalar = length(obj); % number of Scalars to subdivide
scalarIsColumn = iscolumn(obj); % if Scalar is a column, then transpose the subScalars before returning
nSubScalar = size(parmRange,1); % number of subdivisions for each Scalar
% Reparameterize each Scalar into subScalars
switch scalarDimension
    case 1 % Scalars are parameterized arcs
        N = obj.Truncation;
        V = pascal(N); % initialize Pascal matrix of size N
        for iSubScalar = 1:nSubScalar % loop over each row in parmRange first so the computations can be reused
            iSubInterval = parmRange(iSubScalar, :);
            
            if isValid
                iCenter = intval(mean(iSubInterval)); % midpoint
                iRescaling = intval(iCenter - iSubInterval(1)); % radius
                iShiftOperator = intval(zeros(N));
            else
                iCenter = mean(iSubInterval); % midpoint
                iRescaling = iCenter - iSubInterval(1); % radius
                iShiftOperator = zeros(N);
            end
            iCenterPower = iCenter.^(0:N-1);
            for jRow = 1:N
                iShiftOperator(jRow, jRow:end) = V(jRow, 1:N-jRow+1).*iCenterPower(1:N-jRow+1)*iRescaling^(jRow-1);
            end
            
            % loop over each Scalar
            for jScalar = 1:nScalar
                jCoefficient = obj(jScalar).Coefficient'; % pick off coefficients and make into a column
                try
                    subScalarCoefficient = (iShiftOperator*jCoefficient)'; % Coefficients should always be columns
                catch
                    subScalarCoefficient = (iShiftOperator*jCoefficient')'; % but mistakes happen
                end
                subScalar(iSubScalar, jScalar) = Scalar(subScalarCoefficient, scalarBasis);
                
            end
        end
        
    otherwise % spaceParm is a higher dimension surface
        error('not yet implemented for higher dimensions')
end

if scalarIsColumn
    subScalar = subScalar'; % transpose to match layout of Scalars
end

switch nargout
    case 1
        varargout{1} = subScalar; % no coordinate transform
    case 2
        varargout{1} = newMTC; % new material coordinates
        varargout{2} = subScalar;
end
end % end subdivide

% Revision History:
%{
14-Jun-2017 - implemented rigorous subdivision for Taylor parameterizations
13-Aug-2018 - updated for Scalar class
08-Mar-2019 - Removed optional argument for evaluating in time before subdividing. This task should be handled by the Atlas or Chart class.
    Change the variable argument "domain" to be a required input eliminate confusion about what domain we are subdividing with respect to. Optimized
    by keeping the reparameterization information for vectorized calls. Changed the output type to always return a Scalar.
%}

