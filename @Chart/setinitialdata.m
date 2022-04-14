function setinitialdata(obj, initialData, basis)
%SETINITIALDATA - Write initial data to chart as a Scalar
%
%   Inputs:
%       chart - instance of Chart class
%       initialData - (1+d)-dimensional array. (j,:,...,:) is the j^th coordinate in phase space defined as a d-dimsional coefficient array.
%
%   Outputs:
%       chart - chart instance with initial data appended
%
%   Other m-files required: @Chart, @Scalar
%   Subfunctions: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 01-Jul-2018; Last revision: 20-Aug-2018

% append initial data to chart
if isa(initialData, 'double') || isa(initialData, 'intval')
    S = size(initialData);
    %     if isequal(S(1), obj.Dimension(2)) % initialData(j,:) is a phase space coordinate
    coefDataSubs = mat2cell(S(2:end), 1, length(S)-1);
    coefData = mat2cell(initialData, ones(1, S(1)), coefDataSubs{:}); % convert to length-n cell vector
    obj.InitialData = Scalar(coefData, basis);
    %     else
    %         error('The first dimension of initialData must be the same as the phase space')
    %     end
    
elseif isa(initialData, 'Scalar')
    
elseif isa(initialData, 'cell')
    
else
    error('Initial data type not supported')
end
end

% Revision History:
%{
14-Aug-2018 - updated for Scalar class
20-Aug-2018 - updated for changed to Scalar class constructor
19-Mar-2018 - removed error for appending initial data of wrong size. This will make it compatible with automatic differentiation.
%}
