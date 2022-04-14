function copyObj = deepcopy(obj)
%DEEPCOPY - return a deepcopy of a Scalar object
%
%   DEEPCOPY() - The Scalar class is inherited from handle so reassignment of chart objects between variables still share the same instance of the object.
%       For example the code:
%           a = Scalar(with_some_property)
%           b = a
%           b.property = some_different_value
%       results in the corresponding property changing for both a and b since each variable is simply a pointer to the same Scalar object.
%       This function returns a deepcopy which can be changed without altering the original.
%
%   Syntax:
%       copyObj = DEEPCOPY(obj)
%
%   Inputs:
%       obj - A subclass of the Scalar superclass
%
%   Outputs:
%       copyObj - A Scalar with the same subclass as obj and identical properties which references a distinct handle
%
%   Subfunctions: none
%   Classes required: Chart, Scalar
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Mar-2019; Last revision: 09-Mar-2019

warning off MATLAB:structOnObject
if length(obj) == 1 % deepcopy of a single Scalar
    copyObj = Scalar(); % Initialize a new chart of the same type
    objProperties = fields(struct(obj)); % extract all properties of obj including hidden properties
    for iProperty = 1:length(objProperties) % loop over all properties of obj and assign them to copyObj
        thisProperty = objProperties{iProperty};
        copyObj.(thisProperty) = obj.(thisProperty); % assign this property to copyObj
    end
    
else % deepcopy of a vector of Scalars
    nScalar = length(obj);
    copyObj(nScalar) = obj(nScalar).deepcopy(); % initialize copy vector and copy last Scalar
    for iScalar = 1:nScalar-1
        copyObj(iScalar) = obj(iScalar).deepcopy();
    end
    copyObj = reshape(copyObj,size(obj)); % change layout to match input
end
end % end deepcopy

% Revision History:
%{

%}
