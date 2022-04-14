function varargout = plotdecay(obj, varargin)
%PLOTDECAY - plot the coefficient decay of a Scalar
%
%   PLOTDECAY() - A more detailed description of the function
%
%   Syntax:
%       output = PLOTDECAY(input1, input2)
%       [output1, output2] = PLOTDECAY(input1, input2, input3)
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
%   Date: 21-Aug-2018; Last revision: 21-Aug-2018

if nargin == 1
    if obj.Dimension < 3
        if obj.Dimension == 1
        else
            Z = log10(abs(mid(obj.Coefficient)));
            minLevel = max(min(Z(:)), log10(eps)); % minimal value larger than eps
            maxLevel = max(Z(Z(:) < Inf)); % maximal finite values
            normalizeIdx = -Inf < Z & Z < minLevel;
            Z(normalizeIdx) = minLevel; % normalize values smaller than eps
            decRange = maxLevel - minLevel;
            if decRange < 50
                cLevel = minLevel:1:maxLevel;
            else
                cLevel = linspace(minLevel, maxLevel, 50);
            end
            [Y, X] = ndgrid(1:obj.Truncation(1), 1:obj.Truncation(2));
            decayPlot = contourf(X, Y, Z, cLevel);
            colorbar
        end
    else
        error('Specify dimensions to plot the norm')
        
    end
    
else
    
    
end

if nargout > 0
    varargout{1} = decayPlot;
end
end % end plotdecay

% Revision History:
%{
21-Aug-2018 - updated for Scalar class
%}
