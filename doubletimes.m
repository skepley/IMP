function convCoefficient = doubletimes(leftFactor, rightFactor, varargin)
%DOUBLETIMES - Scalar multiplication (linear convolution) for intval coefficients
%
%   DOUBLETIMES() - A more detailed description of the function
%
%   Syntax:
%       output = DOUBLETIMES(input1, input2)
%       [output1, output2] = DOUBLETIMES(input1, input2, input3)
%
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 22-Aug-2018; Last revision: 22-Aug-2018

% unpack truncation settings
if isscalar(leftFactor) || isscalar(rightFactor) % one factor is just a constant
    truncMode = 'Constant';
    convCoefficient = leftFactor*rightFactor;
elseif nargin == 2 || strcmp(varargin{1}, 'Full') % default is to return the full convolution including higher order terms
    truncMode = 'Full';
    dims = size(leftFactor);
    truncation = dims(dims > 1);
    dimension = sum(dims > 1);
else
    truncMode = 'Partial'; % special truncation modes
    fullConvSize = size(leftFactor) + size(rightFactor) - 1;
    [truncSize{1:nargin-2}] = varargin{:}; % A cell array of the form {1:N1, 1:N2,..., 1:Nd} specifying which coefficients to keep
    if isequal(length(truncSize),1) && isequal(length(fullConvSize),2) % Handle MATLAB reporting vectors as dimension 2
        validConv = {truncSize{1} <= max(fullConvSize)};
    else % compare specified truncation indices against size of full convolution in each dimension
        validConv = arrayfun(@(idx)truncSize{idx}(truncSize{idx} <= fullConvSize(idx)), 1:length(fullConvSize), 'UniformOutput', false);
    end
    dimension = length(truncSize);
    
    if ~isequal(size(leftFactor), size(rightFactor)) % truncated convolution requires same size array
        [leftFactor, rightFactor] = Scalar.commonsize(leftFactor, rightFactor);
    end
    truncation = size(leftFactor);
end

% form convolution product
switch truncMode
    case 'Full'
        convCoefficient = convn(leftFactor, rightFactor);
    case 'Partial' % cell array of indices
        switch dimension
            case 0 % multiplication of constants
                convCoefficient = leftFactor*rightFactor;
                
            case 1 % 1-d convolution
                I = truncSize{1};
                leftPad = [zeros(1, length(leftFactor)-min(I)), reshape(leftFactor, 1, [])];
                padConv = conv(leftPad, rightFactor, 'valid');
                convCoefficient = padConv(I-min(I)+1); % row vector of coefficients of discrete convolution
                
                % transpose to match shape of input vectors
                if ~isrow(leftFactor)
                    convCoefficient = convCoefficient.';
                end
                
                
            otherwise % n-d convolution
                % pad enough to make required indices valid and call faster matlab routine
                padLowDim = truncation - cellfun(@(idx)min(idx), validConv);
                padHighDim = cellfun(@(idx)max(idx), validConv) - truncation;
                fullPadSize = truncation + padLowDim + padHighDim;
                insertSubs = Scalar.trunc2subs(padLowDim + 1, padLowDim + truncation); % indices to insert left Factor
                leftPad = zeros(fullPadSize);
                leftPad(insertSubs{:}) = leftFactor;
                padConv = convn(leftPad, rightFactor, 'valid');
                
                % pick out valid indices
                validSubs = cellfun(@(validIdx)validIdx + 1 -min(validIdx), validConv, 'UniformOutput',false);
                convCoefficient = padConv(validSubs{:});
        end
end % switch
end % doubletimes



% Revision History:
%{
22-Aug-2018 - support for arbitrary dimensions
    removed from mtimes method. Functionality now takes double array as input and output. This function is called inside SCalar/mtimes
    but may be called outside as well.
%}
