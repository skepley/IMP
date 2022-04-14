function varargout = local2global(obj, domainPt, varargin)
%LOCAL2GLOBAL - Converts between chart (material) time and real time.
%
%   LOCAL2GLOBAL() - evaluates the linear maps: tau(t):[0,1] ---> [t0, tf] which maps local chart time 
%           coordinates to global time coordinates for the flow and sigma(s):[-1, 1] ---> [a,b] subset [-1, 1] which
%           maps local chart space coordinates to global space coordinates for the manifold.
%
%   Syntax:
%       tau = LOCAL2GLOBAL(obj, t) 
%       t = LOCAL2GLOBAL(obj, tau, -1) evaluates the inverse linear map t(tau): [t0, tf] ---> [-1,1]
%
%   Inputs:
%       domainPt - An ordered tuple of the form (s1, s2,..., s_{d-1}, t)
%       domainIdx - An indexing subset of {1,...,d} describing which coordinates are specified by domainPt. Default is to assume a full vector of space/time coordinates

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Mar-2019; Last revision: 13-Aug-2020

if nargin == 3
    domainIdx = varargin{1};
else
    domainIdx = 1:length(domainPt);
end 

if isequal(domainIdx, -1)
    error('Inverse mapping was removed from this method. This needs to be replaced by a global2local method.')
    %{
OLD CODE
tau = (2*domainPt(end) - sum(obj.TimeSpan))./obj.Tau; % linear interpolant for data [t0,tf], [-1,1] evaluated at t.
    %}
end

% map local (chart) coordinates to global (real) time coordinates
d = obj.Dimension(1); % get dimension of Chart
tau = []; % initialize as empty vectors
sigma = [];

if any(domainIdx < d) % some domain indices are spatial variables
    sigma = .5*(sum(obj.SpatialSpan) + domainPt(1)*diff(obj.SpatialSpan));
end

if ismember(d, domainIdx) % time variable is included in domainIdx
    % OLD CODE - SEE 13-AUG-2020 BUG FIX
    % tau = 0.5*(obj.TimeSpan(2)*(domainPt(end)+1) - obj.TimeSpan(1)*(domainPt(end)-1)); % linear interpolant for data [-1,1], [t0,tf] evaluated at t.
    tau = domainPt(end) * diff(obj.TimeSpan) + obj.TimeSpan(1);  % new version is the map: t --> t * (tf - t0) + t0
   
end
globalCoordinate = [sigma, tau]; % filter and return correct variable indices

switch nargout
    case 0 % just print to screen
        disp(globalCoordinate);
        
    case 1 % output is a d-vector of global [space, time] coordinates
        varargout{1} = globalCoordinate;
        
    case length(globalCoordinate) % unpack global coordinates
        varargout = mat2cell(globalCoordinate, 1, length(globalCoordinate));
        
    otherwise
        error('Number of output arguments must be 0, 1, or d')
end


% Revision History:
%{
08-Mar-2019 - This file replaces the previous version mtCoords.
17-Jun-2019 - Updated documentation. Fixed output for single varargout calls. Removed inverse mapping function. Added support for specifying
    partial indices to return global domain values.
13-Aug-2020 - Fixed a bug in which the local time interval was 
%}
