%UNTITLED3 - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       UNTITLED3 description
%
%   Output:
%       UNTITLED3 output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 04-Aug-2018


%% ================================================== TEST FOURIER/TAYLOR MIXED CONVOLUTION ==================================================


N = 100;

a = intval(randi(10,N,1));
b = intval(randi(10,N,1));
c1 = cconv(mid(a),mid(b),N);
c2 = fouriertimes(a,b);
max(abs(c1-c2))


return

%% ================================================== SECTION 1 ==================================================
a=randi(10,N,1);
b=randi(10,N,1);
for j = 1:3
    cconv(a,b,N);
end

tic
for j=1:1000
    a=randi(10,N,1);
    b=randi(10,N,1);
    cconv(a,b,N);
end
toc

a=randi(10,N,1);
b=randi(10,N,1);
for j = 1:3
    fouriertimes(a,b);
end

tic
for j=1:1000
    a=randi(10,N,1);
    b=randi(10,N,1);
    fouriertimes(a,b);
end
toc



