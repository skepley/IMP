function varargout = eval(obj, evalData, varargin)
%EVAL - Evaluate a obj on space/time domain
%
%   Inputs:
%       obj - instance of Chart class
%       data - m-by-(k+1) of the form [S1,S2,...,Sk,T];
%
%   Outputs:
%       imageData - cell Array of evaluations in each coordinate
%       [Gamma_d,...,Gamma_d] - Evaluations for each scalar coordinate
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jul-2018; Last revision: 10-Jan-2021

%%
% parse input and varargin
p = inputParser();
p.addRequired('obj')
p.addRequired('data')
p.addParameter('GlobalTime', false)
p.addParameter('GlobalSpace', false)
p.parse(obj, evalData, varargin{:})
globalTime = p.Results.GlobalTime;
globalSpace = p.Results.GlobalSpace;

% convert global to local coordinates if needed
if globalTime && globalSpace
    intersectionData = obj.intersectdomain(evalData);
    localTime = (intersectionData(:,end) - obj.TimeSpan(1))./obj.Tau; % convert global time to local time.
    localSpace = (2 * intersectionData(:,1:end-1) - sum(obj.SpatialSpan))./diff(obj.SpatialSpan); % convert global time to local time.
    data = [localSpace, localTime];
    
elseif globalTime
    warning('Chart:eval', 'this may omit endpoints due to floating point error')
    localTime = (evalData(:,end) - obj.TimeSpan(1))./obj.Tau; % convert global time to local time.
    validIdx = (0 <= localTime) & (localTime <= 1); % check data for evaluations which lie in this chart
    data = [evalData(validIdx, 1:end-1), localTime(validIdx)]; % filter out valid evaluations
    
elseif globalSpace
    warning('Chart:eval', 'this may omit endpoints due to floating point error')
    localSpace = (2 * evalData(:,1:end-1) - sum(obj.SpatialSpan))./diff(obj.SpatialSpan); % convert global time to local time.
    validIdx = (-1 <= localSpace) & (localSpace <= 1); % check data for evaluations which lie in this chart
    data = [localSpace(validIdx,:), evalData(validIdx,end)]; % filter out valid evaluations
    
else
    data = evalData; % evalData is local in both space and time
end

%% evaluate data in local coordinates
if ~isempty(data)
    if nargout == obj.Dimension(2)
        for j = 1:obj.Dimension(2) % loop over phase space variables
            varargout{j} = real(obj.Coordinate(j).eval(data));
        end
    elseif nargout == obj.Dimension(2)+1
        for j = 1:obj.Dimension(2) % loop over phase space variables
            varargout{j} = real(obj.Coordinate(j).eval(data));
        end
        varargout{obj.Dimension(2)+1} = data;
    elseif nargout <= 1
        varargout{1} = {};
        for j = 1:obj.Dimension(2)
            varargout{1}{j} = real(obj.Coordinate(j).eval(data));
        end
    else
        error('nargout for eval all should be 1 or d')
    end
    
else
    if nargout == 1
        varargout{1} = {};
    else
        
        for j = 1:nargout
            varargout{j} = [];
        end
    end
end %  if
end %  eval

% Revision History:
%{
23-Mar-2019 - Added evaluation in global space/time coordinates
10-Jan-2021 - Renamed the input variable to "evalData" while the global
space/time variable keeps the name "data". This makes it easier to debug
and troubleshoot methods related to evaluation. 
%}
