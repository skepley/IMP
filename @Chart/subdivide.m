function subChart = subdivide(obj, parmRange)
%SUBDIVIDE - Subdivide a chart into smaller charts which are piecewise parameterizations of the same manifold
%
%   SUBDIVIDE() - A more detailed description of the function
%
%   Syntax:
%       output = SUBDIVIDE(input1, input2)
%       [output1, output2] = SUBDIVIDE(input1, input2, input3)
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
%   Date: 08-Mar-2019; Last revision: 08-Mar-2019
%
% TODO: Take advantage of improved error bounds for subcharts in the interior of the spatial domain

if isequal(parmRange, obj.SpatialSpan) % no subdivision is required
    subChart = obj;
else
    nSubChart = size(parmRange,1);
    coordinateIsColumn = iscolumn(obj.Coordinate); % check layout of Coordinate vector
    [subDomain, subCoordinate] = obj.Coordinate.subdivide(parmRange, obj.SpatialSpan); % Subdivision of each Chart coordinate
    for iSubChart = 1:nSubChart
        subChart(iSubChart) = obj.deepcopy(); % initialize each subchart
        subChart(iSubChart).SpatialSpan = subDomain(iSubChart,:); % update the subdomain
        if coordinateIsColumn
            subChart(iSubChart).Coordinate = subCoordinate(:, iSubChart);
        else
            subChart(iSubChart).Coordinate = subCoordinate(iSubChart, :);
        end
    end
end
end % end subdivide

% Revision History:
%{

%}
