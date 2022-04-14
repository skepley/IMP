function varargout = solveradiipoly(coef, varargin)
%SOLVERADIIPOLY - a rigorous bound on smallest root of a quadratic polynomial specified (usually) as [Z2, -(1-Z1), Y0].
%
%   r1 = SOLVERADIIPOLY(coef) returns a rigrous enclosure on the smallest root for the polynomial.
%
%   r1 = SOLVERADIIPOLY(coef, true) prints accompanying midrad values for the root and compares to matlabs built in solver.
%
%   [r1, r2]  = SOLVERADIIPOLY(coef, true) returns rigorous enclosure for both roots
%
%   Inputs:
%       coef - intval
%       printStats - logical
%
%   Outputs:
%       r1 - intval
%       r2 - intval
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: INTLAB
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Aug-2018; Last revision: 19-Aug-2018

if nargin >1
    printStats = varargin{1};
else
    printStats = false;
end

polyCoef = intval(coef);
radiiPoly = polynom(polyCoef, 'r'); % radii polynomial
try
    rGuess = roots(mid(polyCoef)); % numerical estimate for true roots.
catch ME
    disp('here')
end
if nargout ==1 % return only smallest root
    rTrue = verifypoly(radiiPoly, min(rGuess));
    varargout{1} = rTrue;
    
else % return both roots
    rTrue(1) = verifypoly(radiiPoly, min(rGuess));
    rTrue(2) = verifypoly(radiiPoly, max(rGuess));
    varargout{1} = rTrue(1);
    varargout{2} = rTrue(2);
end

if printStats % print statistics
    fprintf('\n')
    fprintf('Rpoly midpoint: %.3e, %.3e, %.3e \n',mid(polyCoef))
    fprintf('Rpoly radius: %.3e, %.3e, %.3e \n',rad(polyCoef))
    fprintf('intval root: %.3e, %.3e \n',mid(rTrue(1)),rad(rTrue(1)))
    if nargout ==1
        fprintf('Matlab root solver: %.3e \n', min(rGuess))
    else
        fprintf('Matlab root solver: %.3e, %.3e \n', rGuess)
    end
    fprintf('\n')
end
end

% Revision History:
%{
19-Aug-2018 - Added documentation, updated for Scalar class and changed output to intval.
%}
