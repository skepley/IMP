function varargout = bestfitdecay(obj)
%BESTFITDECAY - Returns the best approximate geometric sequence fitted to the given coefficents by linear least squares regression.
%
%   r = BESTFITDECAY(obj) is a vector of ratios, [r_1,...,r_d] with length obj.Dimension. These ratios are the best fit of Scalar.Coefficient to
%       a model geometric series with geometric rate r_i in the ith direction obtained by log-linear regression i.e. this is the slope of the best
%       fit line for the data given by log(Scalar.Coefficient)
%
%   [r, A] = BESTFITDECAY(obj) is the best fit ratios and the best fit initial term i.e. this is both the slope and the intercept for the line of best
%       fit. 
%
%       1D - If the coefs are of the form a = (A, Ar, Ar^2, ... Ar^{n-1}), then log a lies on the line y(x) = log A + log(r)*x.
%           The best fit (log) linear line is the minimizer of the sum of squared residuals S = sum(r_i^2) = sum(log(a_i) - beta_0 - i*beta_1)^2.
%       2D - The returned ratios minimize S = sum_j sum_i (log(a_ij) - beta_0 - j*beta_1 - i*beta_2)^2 where j is the column index, and i is the row index.
%
%   Inputs:
%       obj - Scalar of dimension d
%
%   Outputs:
%       geometricDecay - varies
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jan-2018; Last revision: 1-Feb-2019

% TODO: Replace 1 dimensional code with the MATLAB builtin polyfit

if length(obj) > 1 % vectorize method
    for j = 1:length(obj)
        if isequal(nargout,1)
            varargout{1}(j,:) = bestfitdecay(obj(j)); % output slope
        else
            [varargout{1}(j,:), varargout{2}(j,:)] = bestfitdecay(obj(j)); % output slope and intercept
        end
    end
else
    
    if ~strcmp(obj.Basis, 'Taylor')
        warning('bestfitdecay - fitting to a geometric series is only reasonable for Taylor series')
    end
    
    if obj.Dimension == 1 % coef is a 1D coefficient vector
        logCoefficient = log(abs(obj.Coefficient)); % best fit geometric sequence is linear best fit for log(sequence)
        nonZeroIdx = logCoefficient > -Inf;
        Y = logCoefficient(nonZeroIdx); % remove zero coefficients from regression data
        N = length(logCoefficient);
        X = reshape(0:N-1,size(logCoefficient));
        X = X(nonZeroIdx); % remove zero coefficients from regression nodes
        Xsym = [N,sum(X);sum(X),sum(X.^2)]; % get symmetric part of X
        XsymY = [sum(Y);sum(X.*Y)];
        beta = Xsym\XsymY;
        
    elseif obj.Dimension == 2
        logCoefficient = reshape(log(abs(obj.Coefficient))', [], 1);
        nonZeroIdx = logCoefficient > -Inf;
        Y = logCoefficient(nonZeroIdx);
        [M,N] = size(obj.Coefficient);
        S = 0:N-1;
        T = (0:M-1);
        
        % 2D least squares regression
        X = zeros(M*N,3);
        X(:,1) = ones(M*N,1);
        X(:,2) = repmat(S',M,1);
        X(:,3) = reshape(repmat(T,N,1),[],1);
        X = X(nonZeroIdx,:);
        XTX = X'*X; % get symmetric part of X
        
        XTY = X'*Y;
        beta = XTX\XTY;
    else
        error('not implemented for higher dimensions')
    end
    
    % output
    varargout{1} = exp(beta(2:end)); % output slopes via best fit geometric ratio
    if isequal(nargout, 2)
        varargout{2} = exp(beta(1)); % constant multiple for geometric sequence
    end
    
end % bestfitdecay

% Revision History:
%{
13-Aug-2018 - updated for Scalar class
29-Aug-2018 - fixed log-linear regression to handle the case when some coefficients are zero
16-Jan-2019 - fixed a bug in the 1-D version when some coefficients were zero.
1-Feb-2019 - Added support for getting the intercept of the best fit line. 
%}

