function taylor_coefs = bivar_polycomp(P,theta1,theta2)
% Computes P(s) given coefficients for P(theta1,theta2) and theta1(s), theta2(s). 
m = size(P,1);
n = size(P,2);
T1 = polycomp(theta1,m-1);
T2 = polycomp(theta2,n-1);

if isa(T1,'cell') % theta has intval coefficients 
    deg = length(T1) + length(T2) - 2;
    taylor_coefs = intval(zeros(1,deg+1));
    for j=1:m
        for k = 1:n
            addPoly = P(j,k)*T1{k}*T2{j};
            addto = flip(addPoly.c);
            taylor_coefs(1:length(addto)) = taylor_coefs(1:length(addto)) + addto;
        end
    end
    
else   % theta has double coefficients
    deg = size(T1,2) + size(T2,2) - 2;
    taylor_coefs = zeros(1,deg+1);   
    if isa(P,'intval')
        taylor_coefs = intval(taylor_coefs);
    end    
    for j=1:m
        for k = 1:n
            addto = P(j,k)*conv(T1(k,:),T2(j,:));
            taylor_coefs(1:length(addto)) = taylor_coefs(1:length(addto)) + addto;
        end
    end
end

end