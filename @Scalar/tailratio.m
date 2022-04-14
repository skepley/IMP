function tailRatio = tailratio(obj, tailCutoff)
%TAILRATIO - computes the norm of the higher order terms of a Scalar relative to Scalar's total norm
%
%   objTailRatio = TAILRATIO(obj, N) returns a value in [0,1]. N = [N1,...,Nd] is an integer vector of cutoffs which determine the "tail".
%       Coefficient indices which exceed N are members of the tail. Their norm is divided by the norm of all coefficients.
%
%   Inputs:
%       obj - Scalar
%       tailCutoff - Integer vector of length d
%
%   Outputs:
%       tailRatio - double
%
%   Subfunctions: NORM
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 23-Apr-2017; Last revision: 13-Aug-2018

if ~strcmp(obj(1).Basis, 'Taylor')
    warning('tailratio - may not be appropriate for non-Taylor basis')
end

objDimsidx = Scalar.dims(obj(1).Coefficient);
if objDimsidx > 1
    error('This function has not been updated for higher dimensional data')
end
Sz = size(obj(1).Coefficient);
Sz = Sz(objDimsidx);
if any(tailCutoff < 1) % Non-tail terms specified as fraction of total terms.
    fractionTerms = tailCutoff < 1;
    tailCutoff(fractionTerms) = round(tailCutoff(fractionTerms).*Sz(fractionTerms));
    tailCutoff(tailCutoff < 1) = ones(size(tailCutoff < 1));
end

if length(obj) > 1 % vectors of Scalars
%     getRatio = @(j)obj(j).tailratio(tailCutoff);
    tailRatio(length(obj)) = obj(end).tailratio(tailCutoff);
    for j = 1:length(obj)-1
       tailRatio(j) = obj(j).tailratio(tailCutoff);
    end
    tailRatio = reshape(tailRatio, size(obj));
else
    switch length(tailCutoff)
        case 1 % obj is 1d
            objTail = Scalar(obj.Coefficient(tailCutoff:end), obj.Basis);

        case 2 % obj is 2d
            tailCoefficient = obj.Coefficient;
            tailCoefficient(1:tailCutoff(1)-1, 1:tailCutoff(2)-1) = zeros(1:tailCutoff(1)-1, 1:tailCutoff(2)-1);
            objTail = Scalar(tailCoefficient, obj.Basis);

        case 3 % obj is 3d
            objTail = Scalar(obj.Coefficient(tailCutoff(1):end, obj.Basis, tailCutoff(2):end, tailCutoff(3):end));

        otherwise
            error('tailRatio not implemented for Scalar with dimension greater than 3')
    end
    tailNorm = objTail.norm();
    tailRatio = tailNorm./(obj.norm - abs(obj.Coefficient(1))); % 9 July 2019 - Subtract off the constant term  
end
if isa(tailRatio,'intval')
    tailRatio = sup(tailRatio);
end
end % tailratio

% Revision History:
%{
19-Aug-2017 - support for intval coefficients
13-Aug-2018 - updated for Scalar class
21-Mar-2018 - fixed bugs from previous update
09-Jul-2019 - Changed the computed ratio to exclude the constant term since this can hide a "fat" tail if it is small relative to the constant. 
%}
