function varargout = truncatepoly(poly, varargin)
%TRUNCATEPOLY - returns the reductum (truncation) of an intlab polynomial
%
%   pR = TRUNCATEPOLY(poly, [1,M]) returns the polynomnial with terms indexed by varIdx a univariate intlab polynom.
%
%   p_M = TRUNCATEPOLY(poly, M) returns the M^th coefficient of the polynomial polynomnial with terms indexed by varIdx a univariate intlab polynom.
%
%   pT = TRUNCATEPOLY(poly, var1Idx, var2Idx) truncates an intlab polynom in 2 variables.
%
%   [pT, T] = TRUNCATEPOLY(p, varIdx) returns a rigorous enclosure of the ell_1 norm of the discarded coefficients of p.
%
%   Input:
%   poly - IntLab polynom class object
%   M - single integer
%   [m, M] - length 2 integer
%
%   Output:
%   pR - IntLab polynom class object
%   p_M - 1-by-1 intval
%   T - 1-by-1 double
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: IntLab
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jan-2018; Last revision: 18-Jul-2018


truncVar = true(1,length(poly.v));
if ~isequal(length(truncVar),nargin -1)
    error('Must pass indices for each variable')
end

switch nargin
    case 2 % univariate polynomial
        idx = varargin{1}; % index of the form M or [1,M] or [m,M]
        idx = idx(idx <= 1 + poly.e); % remove indices greater than length of coefficients
        if isempty(idx)
            truncatedPoly = poly;
        else
            flipCoef = flip(poly.c); % coefficients in ascending order
            truncCoef = flipCoef(idx); % pick off indices
            truncatedPoly = polynom(truncCoef,idx-1,poly.v(truncVar));  % return truncated polynomial
            if nargout ==2 % return bounds on tail terms
                tailIdx = true(1,1 +poly.e);
                tailIdx(idx) = false;
                tailCoef = flipCoef(tailIdx);
                varargout{2} = sum(abs(tailCoef));
            end
        end
        
    case 3 % multi-variate polynomial
        %         idx1 = varargin{1}; % index of the form N or [1,N] or [n,N]
        %         idx2 = varargin{2}; % index of the form M or [1,M] or [m,M]
        
        coefIdx = true(size(poly.e,1),1);
        for j = 1:nargin -1
            idx = varargin{j};
            if ~isempty(idx) % truncate nothing if idx = []
                if length(idx) ==1 % idx = N
                    varIdx = poly.e(:,j) == idx-1;
                    truncVar(j) = false; %truncatedPoly is constant with respect to this variable.
                elseif isequal(idx(1),1) % index type [1,N]
                    varIdx = poly.e(:,j) < idx(2);
                else
                    varIdx = poly.e(:,j) < idx(2) & idx(1)-1 <= poly.e(:,j);
                end
                coefIdx = coefIdx & varIdx;
            end
        end
        truncatedPoly = polynom(poly.c(coefIdx),poly.e(coefIdx,truncVar),poly.v(truncVar));
end % switch
varargout{1} = truncatedPoly;
end % truncatepoly

% Revision History:
%{
08-Aug-2018 - updated for Scalar class
19-Aug-2018 - fixed constant projection case, added documentation
%}



