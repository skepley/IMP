function powerArray = polycomp(p,m)
% p is a polynomial of degree d (coefs in INCREASING order), returns a matrix whose (k+1)^th row is p(s)^k in INCREASING order. 
deg = length(p) - 1;
if isa(p,'intval')
    P = polynom(fliplr(p)); % intLab has coefficients in decreasing order
    powerArray = cell(m,1);
    powerArray{1} = polynom(1);
    for j = 2:m+1
        powerArray{j} = P*powerArray{j-1};
    end
else
    M = zeros(m,deg*m + 1);
    M(1,1:deg+1) = p;
    for j = 2:m
        row = conv(p,M(j-1,:));
        M(j,:) = row(1:size(M,2));
    end
    powerArray = [zeros(1,size(M,2));M];
    powerArray(1,1) = 1;
end
end