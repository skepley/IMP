function copyObj = deepcopy(obj)
%DEEPCOPY - return a deepcopy of a Chart object
%
%   DEEPCOPY() - The Chart class is inherited from handle so reassignment of chart objects between variables still share the same instance of the object.
%       For example the code:
%           a = Chart(with_some_property)
%           b = a
%           b.property = some_different_value
%       results in the corresponding property changing for both a and b since each variable is simply a pointer to the same Chart object.
%       This function returns a deepcopy which can be changed without altering the original.
%
%   Syntax:
%       copyObj = DEEPCOPY(obj)
%
%   Inputs:
%       obj - A subclass of the Chart superclass
%
%   Outputs:
%       copyObj - A Chart with the same subclass as obj and identical properties which references a distinct handle
%
%   Subfunctions: none
%   Classes required: Scalar, Chart
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Mar-2019; Last revision: 09-Mar-2019

warning off MATLAB:structOnObject
if length(obj) == 1 % deepcopy of a single Chart
    copyObj = obj.ClassConstructor(); % Initialize a new chart of the same type
    objProperties = fields(struct(obj)); % extract all properties of obj including hidden properties
    for iProperty = 1:length(objProperties) % loop over all properties of obj and assign them to copyObj
        thisProperty = objProperties{iProperty};
        copyObj.(thisProperty) = obj.(thisProperty); % assign this property to copyObj
    end
    % Scalars are also inherited from handle so the coordinates must be deepcopied from the Scalar method
    copyObj.Coordinate = deepcopy(obj.Coordinate);
    
else % deepcopy of a vector of Charts
    nChart = length(obj);
    copyObj(nChart) = obj(nChart).deepcopy(); % initialize copy vector and copy last Chart
    for iChart = 1:nChart-1
        copyObj(iChart) = obj(iChart).deepcopy();
    end
    copyObj = reshape(copyObj,size(obj)); % reshape to same layout as input Chart
end
end % end deepcopy

% Revision History:
%{

%}
