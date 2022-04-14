function fout = intlabpoly2scalar(poly, basis, varargin)
%INTLABPOLY2BASCALAR - Convert an intlab polynomial class to a BAscalar class 
%
%   Syntax:
%       output = INTLABPOLY2BASCALAR(input1, input2)
%       output = INTLABPOLY2BASCALAR(input1, input2, input3)
%
%   Description:
%       INTLABPOLY2BASCALAR() - description
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
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Jul-2018; Last revision: 19-Jul-2018

%% 
dim = length(poly.v);
switch dim
    case 1
        fout = Scalar(flip(poly.c), basis); % reverse order to put into ascending order
    otherwise
        error('intlabpoly2scalar not yet implemented for higher dimensions')
end

end % end intlabpoly2BAscalar

