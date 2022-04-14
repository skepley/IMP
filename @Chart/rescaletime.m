function rescaletime(obj,varargin)
%RESCALETIME - Automatic time rescaling for Charts
%
%   RESCALETIME() - A more detailed description of the function
%
%   Syntax:
%       output = RESCALETIME(input1, input2)
%       [output1, output2] = RESCALETIME(input1, input2, input3)
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
%   Date: 22-Mar-2019; Last revision: 22-Mar-2019
if nargin > 1
    lastNorm = varargin{1};
else
    lastNorm = eps(1);
end

rowDecay = obj.Coordinate.decay;
maxDecay = max(rowDecay(end,:));
L = (lastNorm/maxDecay)^(1/(obj.Truncation(1)-1));
scaleBy = mid(L*obj.Tau);
obj.scaletime(scaleBy);
end % end rescaletime

% Revision History:
%{

%}
