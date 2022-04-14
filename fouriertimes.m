function convObj = fouriertimes(leftObj,rightObj)
%FOURIERTIMES - Computes the product (circular convolution) of two truncated discrete Fourier series
%
%   Syntax:
%       output = FOURIERTIMES(input1, input2)
%       output = FOURIERTIMES(input1, input2, input3)
%
%   Description:
%       FOURIERTIMES() - description
%
%   Inputs:
%       leftObj - Fourier coefficients indexed as [a_0,...,a_{N-1}]
%       leftObj - Fourier coefficients indexed as [b_0,...,b_{N-1}]
%
%   Outputs:
%       convObj -  Fourier coefficients indexed as [c_0,...,c_{N-1}] corresponding to a*b
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Aug-2018; Last revision: 04-Aug-2018

%%
dim = sum(size(leftObj) > 1);

switch dim
    case 1
        N = length(leftObj);
        if length(rightObj) ~= N
            warning('not tested for Fourier coefficients of different length')
        end

        if isa(leftObj, 'intval') || isa(rightObj, 'intval') % intval circular convolution
            % build circulant convolution matrix
            idx = [N, 1:N-1];
            indexMatrix = flipud(toeplitz([idx(1) idx(end:-1:2)],idx));
            coefMatrix = reshape(leftObj,[],1)*reshape(rightObj,1,[]);  % outer product
            convObj = full(sparse(indexMatrix,1,coefMatrix));
        else
            convObj = cconv(leftObj,rightObj,N); % Signal processing builtin is FFT based so its much faster than mine but only works for doubles
        end

    otherwise
        error('not implented yet')

end
end %  fouriertimes

