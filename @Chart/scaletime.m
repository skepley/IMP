function scaletime(obj, scaleBy)
%SCALETIME - Rescale time for a given chart based on single-step options
%
%   Description:
%       SCALETIME() - determine time scaling, L such that the time-1 map of the rescaled flow is equivalent to the time-L map of the given flow.
%
%   Inputs:
%       chart - instance of the Chart class
%       scaleBy - double which represents the new time tau (in global time) for the chart. This should be strictly positive unless the rescaling
%             is intended to change the flow direction as well.
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Mar-2017; Last revision: 01-Jul-2018


%%
if ~strcmp(obj.Weight,'ones') % Assumption is that ell^1 weights should always be normalized to 1
    error('Weights other than one are not supported.')
end

if strcmp(obj.NumericalClass, 'intval')
    warning('Interval coefficients should never be rescaled this way.')
end

if ~isequal(sign(scaleBy), sign(obj.Tau))
    warning('Scaling time by a negative value will reverse the time direction of the integration. Be sure this is intentional')
end

obj.Coordinate = obj.Coordinate.scaletime(scaleBy/obj.Tau); % call Scalar/scaletime on phase coordinates
obj.Tau = scaleBy; % update chart timestep
obj.TimeSpan(2) = obj.TimeSpan(1) + obj.Tau; % update chart timespan

end %  scaletime

% Revision History:
%{
5-July-2019 - Fixed a bug which was flipping the sign of the timestep for backward time charts. Chart.Tau is now always equal to (tFinal - tINitial).
    In particular, integrating backward in time is identical to forward with the difference captured by the Tau property.
%}
