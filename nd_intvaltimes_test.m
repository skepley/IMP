%ND_INTVALTIMES_TEST - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       ND_INTVALTIMES_TEST description
%
%   Output:
%       ND_INTVALTIMES_TEST output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 21-Aug-2018; Last revision: 21-Aug-2018
% 
% clc
% a = 1:5;
% I = [2,4];
% m = min(I);
% padA = [zeros(1, length(a)-m), reshape(a, 1, [])];
% b = 5:-1:1;
% C = conv(padA, b, 'valid');
% C(I-min(I)+1)
% c = conv(a,b)


a = randi(5,1,5);
b = randi(5,1,2);
cf = convn(a,b)

convIdx = {1:20,4:5};
convDim = size(a) + size(b) - 1;
validIdx = arrayfun(@(idx)convIdx{idx}(convIdx{idx} <= convDim(idx)), 1:length(convDim), 'UniformOutput', false);
c1 = doubletimes(a,b, convIdx{:})
cf(validIdx{:})


return
% trunc1
i1 = [3,5];
p1 = i1 - 1;
fullSize = size(a) + p1;
% padSubs = arrayfun(@(idx)p1(idx)+1:fullSize(idx), 1:2, 'UniformOutput',false);
clear padA
% padA(padSubs{:}) = a
padA = zeros(fullSize -1);
padA(1:3,1:5) = a;
c1 = convn(padA, b, 'valid')
% cf(1:i1(1),1:i1(2))
cf



% doubletimes(a,b,1:3,1:5)


return
doubletimes(mid(a),mid(b),1:3,1:20)
return

tic
for j = 1:100
    intvaltimes(A,B,'Full');
end
toc

tic 
for j = 1:100
    doubletimes(mid(A.Coefficient), mid(B.Coefficient), 'Full');
end
toc




return

C = intval(Scalar.randi(5,40,20));
tic
for j = 1:25
    c = mtimes(A,B,40,1:20);
    C.append(c);
end
toc

CC = C.Coefficient;
tic 
for j = 1:25
    c = intvaltimes(A.Coefficient, B.Coefficient,40,1:20);
    CC(41,:) = c;
end
toc



