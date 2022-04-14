function [output1,output2] = growatlas(input1,input2,input3,varargin)
%GROWATLAS - Main loop for iteratively generating an atlas by growing the boundary
%
%   GROWATLAS() - Iteratively calls growboundary on the leaves of the chart tree for an atlas and appends the new leaves. This continues until a
%       stopping condition has been reached for every leaf in the tree.
%
%   Syntax:
%       output = GROWATLAS(input1, input2)
%       [output1, output2] = GROWATLAS(input1, input2, input3)
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
%   Date: 07-Mar-2019; Last revision: 07-Mar-2019

%% parse input
p = inputParser;
addRequired(p,input1)
addRequired(p,input2)
addRequired(p,input3)
addParameter(p,'Parameter1',default1)
addParameter(p,'Parameter2',default2)
addParameter(p,'Parameter3',default3)

parse(p,input1,input2,input3,varargin{:})
parameter1 = p.Results.Parameter1;
parameter2 = p.Results.Parameter2;
parameter3 = p.Results.Parameter3;

end % end growatlas

% Revision History:
%{

%}
