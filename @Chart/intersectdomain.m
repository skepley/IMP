function data = intersectdomain(obj, data)
%INTERSECTDOMAIN - returns the subset of evalData which lies in the (global) domain of this chart
%
%   INTERSECTDOMAIN() - Given a chart which parameterizes a function on a rectangle, [t0, t1] X [s0, s1] where [s0, s1] is a subset of [-1,1].
%   evalData is a k-by-2  matrix of points in R^2. This function returns the subset of points which lie in the rectangle.
%
%   Syntax:
%       output = INTERSECTDOMAIN(input1, input2)
%       [output1, output2] = INTERSECTDOMAIN(input1, input2, input3)
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
%   Date: 02-Apr-2019; Last revision: 02-Apr-2019

% filter by time coordinate
globalTime = data(:, end);
t0 = min(obj.TimeSpan);
t1 = max(obj.TimeSpan);
validIdx = (t0 <= globalTime) & (globalTime <= t1); % check data for evaluations which lie in this chart
data = [data(validIdx, 1:end-1), globalTime(validIdx)]; % filter out valid evaluations

% filter by space coordinate
% globalSpace = data(:, 1:end); THIS IS A BUG I THINK. FIXED JUN 18 2019. 
globalSpace = data(:, 1:end-1);
s0 = obj.SpatialSpan(1);
s1 = obj.SpatialSpan(2);
validIdx = (s0 <= globalSpace) & (globalSpace <= s1); % check data for evaluations which lie in this chart
data = [globalSpace(validIdx, :), data(validIdx, end)]; % filter out valid evaluations

end % end intersectdomain

% Revision History:
%{

%}
