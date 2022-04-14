function varargout = ellonebox(obj)
%ELLONEBOX - returns a 1-by-n interval vector which bounds the chart in phase space using coarse ell^1 bounds.

%   Syntax:
%       chart.ELLONEBOX() - Append ell^1 box to chart.ellOneBox property
%       box = chart.ELLONEBOX() - Append and return box
%
%   Description:
%       ELLONEBOX() - description
%
%   Inputs:
%       chart - instance of the Chart class
%
%   Outputs:
%       box - 1-by-n intval
%
%   Other m-files required: @Chart, @Scalar
%   Subfunctions: Scalar/norm
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 01-Jul-2018; Last revision: 01-Jul-2018

if length(obj) > 1 % vectorized version
    for j = 1:length(obj)
        obj(j).ellonebox();
    end
    
else % compute box for a single chart
    center = arrayfun(@(j)obj.Coordinate(j).Coefficient(1,1), (1:obj.Dimension(2))'); % constant term in chart expansion
    fullNorm = obj.Coordinate.norm;
    
    if obj.IsValid
        radius = fullNorm - abs(center) + obj.StepError; % chart.StepError = 0 for nonrigorous timesteps.
    else
        radius = fullNorm - abs(center) + 1e-14; % add a small floating point padding to subdue roundoff error
    end
    box = midrad(center, radius)';
    %     chart.ellOneBox = box;
end

varargout{1} = box; % return box if argout
end

% Revision History:
%{
21 May 2019 - Updated for new Scalar and Chart classes. Added support for nonrigorous chart ellonebox enclosures.
13 Jun 2019 - Fixed bug when computing ell-one boxes for floating point Charts
%}

