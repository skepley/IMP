function convObj = fouriertaylortimes(obj, rightObj)
%FOURIERTAYLORTIMES - mixed convolution for Fourier-Taylor series
%
%   Description:
%       convObj = FOURIERTAYLORTIMES(obj, rightObj) returns the convolution of two Scalars. obj and rightObj must
%       have the same basis type. Convolution is circular in Fourier directions and linear in Taylor directions.
%
%   Inputs:
%       obj -  mixed basis Scalar
%       rightObj - mixed basis Scalar
%
%   Outputs:
%       convObj - Scalar of same basis type
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 13-Aug-2018

%   TODO:
%   add support for truncated convolutions
%   integrate into mtimes method
%   add option for FFT based convolution

switch obj.Dimension
    case 2
        [M,N] = size(obj.Coefficient);
        if ~isequal([M,N], size(rightObj.Coefficient))
            error('mixed basis convolution not tested for non-uniform sized Scalars')
        end

        if strcmp(obj.NumericalClass, 'double')
            linearConv = conv2(obj.Coefficient, rightObj.Coefficient); % full linear convolution in both directions
            circularConv = linearConv(1:M,1:N); % truncate in the Taylor (time) direction
            circularConv(:,1:N-1) = circularConv(:,1:N-1) + linearConv(1:M, N+1:end); % circular convolution in Fourier (spatial) direction
            convObj = Scalar(circularConv, obj.Basis);
        elseif strcmp(obj.NumericalClass, 'intval')

        end
    otherwise
        error('not implemented yet')

end
end %  fouriertaylortimes

% Revision History:
%{
13-Aug-2018 - updated for Scalar class
%}
