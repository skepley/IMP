function taylorCoefs = polycompose(P, theta1, theta2)
%POLYCOMPOSE - compose two polynomials on the level of Taylor coefficients
%
%   POLYCOMPOSE() - A more detailed description of the function
%
%   Syntax:
%       output = POLYCOMPOSE(input1, input2)
%       [output1, output2] = POLYCOMPOSE(input1, input2, input3)
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
%   Date: 04-Apr-2019; Last revision: 04-Apr-2019

% Computes P(s) given coefficients for P(theta1,theta2) and theta1(s), theta2(s).

nDim = nargin-1; % number of variables for P
N = size(P);

switch nDim
    case 2
        T1 = polycomp(theta1,N(1)-1);
        T2 = polycomp(theta2,N(2)-1);
        
        if isa(T1,'cell') % theta has intval coefficients
            deg = length(T1) + length(T2) - 2;
            taylorCoefs = intval(zeros(1,deg+1));
            for j=1:N(1)
                for k = 1:N(2)
                    addPoly = P(j,k)*T1{k}*T2{j};
                    addto = flip(addPoly.c);
                    taylorCoefs(1:length(addto)) = taylorCoefs(1:length(addto)) + addto;
                end
            end
            
        else   % theta has double coefficients
            deg = size(T1,2) + size(T2,2) - 2;
            taylorCoefs = zeros(1,deg+1);
            if isa(P,'intval')
                taylorCoefs = intval(taylorCoefs);
            end
            for j=1:N(1)
                for k = 1:N(2)
                    addto = P(j,k)*conv(T1(k,:),T2(j,:));
                    taylorCoefs(1:length(addto)) = taylorCoefs(1:length(addto)) + addto;
                end
            end
        end
        
    otherwise
        error('not implemented')
end

end % end polycompose

% Revision History:
%{

%}
