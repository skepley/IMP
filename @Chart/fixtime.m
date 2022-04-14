function evalObj = fixtime(obj, t)
%FIXTIME - Evaluates first dimension of Chart and returns a new Chart with lower dimensional domain
%
%   FIXTIME() - Returns a chart of one lower dimension by evaluating the first dimension in every Scalar coordinate at the specified (chart) time. The
%   new chart takes the old chart as its parent.
%
%   Syntax:
%       g = FIXTIME(G, tau) returns a chart of dimension (d-1) obtained by evaluating the d-dimensional chart G at t=tau in the first variable.
%
%   Inputs:
%       obj - Chart object
%       t - A real number in the interval [-1,1] (for Taylor basis)
%
%   Outputs:
%       evalObj - A chart of the same type as obj
%
%   Subfunctions: none
%   Classes required: Chart, Scalar
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Mar-2019; Last revision: 08-Mar-2019

if ~isequal(length(obj.TimeSpan),2)
    error('This chart has not been advected')
end
evalObj = obj.deepcopy(); % Initialize a new chart with same properties as obj
evalObj.Basis = obj.Basis(2:end); % remove first dimensions basis
evalObj.Coordinate = fixtime(obj.Coordinate, t); % evaluate the Scalar coordinates in 1st dimension
evalObj.InitialData = evalObj.Coordinate.deepcopy(); % copy initial coordinates as the initial data Scalar
evalObj.TimeSpan = obj.local2global(t, obj.Dimension(1)); % initial time of new chart is evaluation time of the old chart
evalObj.Dimension(1) = obj.Dimension(1)-1; % drop dimension by one due to evaluation
evalObj.ParentHandle = obj; % define original chart as the parent of the new boundary chart
evalObj.Generation = obj.Generation + 1; % advance to next generation
evalObj.Tau = []; % clear Tau property
evalObj.InitialError = obj.ErrorBound; % pass obj error as initial error to the boundary chart
evalObj.StepError = []; % pass obj error as initial error to the boundary chart
% evalObj.ErrorBound = [];
end % end fixtime

% Revision History:
%{

%}
