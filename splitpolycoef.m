function varargout = splitpolycoef(poly,varargin)
% poly is an intlab polynomial in 1 or 2 variables, p(s) or p(s,t)
% idx is a cell array or vector which describes the indices at which to split the coefficients
% truncated polynomial (i.e. Degree + 1)

% arg types: 
% n: 1-by-1 integer returns the n^th coefficient in that dimension
% [m,n]: 1-by-2 integer returns coefficients between m^th and n^th (inclusive)


%   Subfunctions: see truncatepoly
%   Classes required: none
%   Other m-files required: IntLab
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jan-2018; Last revision: 18-Jul-2018


truncVar = true(1,length(poly.v));
if ~isequal(length(truncVar),nargin -1)
    error('Must pass indices for each variable')
end

if isequal(poly.e, 0)
    varargout{1} = poly.c;
    return
end

switch nargin
    case 2 % univariate polynomial
        idx = varargin{1}; % index of the form M or [1,M] or [m,M]
        if length(idx) ==1 % idx of the form M
            truncatedCoef = poly.c(end +1 - idx); % switch to ascending order
        elseif isequal(idx(1),1) % index type [1,M]
            flipCoef = flip(poly.c); % coefficients in ascending order
            coefIndex = 1:idx(2);
            coefIndex = coefIndex(coefIndex <= 1 + poly.e); % remove indices greater than length of coefficients 
            truncatedCoef = flipCoef(coefIndex); % pick off indices
        end
        
        if nargout ==2 % return bounds on tail terms
            tailIdx = true(1,1 +poly.e);
            tailIdx(coefIndex) = false;
            tailCoef = flipCoef(tailIdx);
            varargout{2} = sum(abs(tailCoef));
        end
       
    case 3 % multi-variate polynomial
%         idx1 = varargin{1}; % index of the form N or [1,N] or [n,N]
%         idx2 = varargin{2}; % index of the form M or [1,M] or [m,M]
        error('This needs to be fixed before using. It was just copied from truncatepoly.m')
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
end
varargout{1} = truncatedCoef;
end
		
		
		

