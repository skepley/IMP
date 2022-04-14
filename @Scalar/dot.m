function dotProduct = dot(obj, rightFactor, varargin)
%DOT - compute the dot product of two Scalars
%
%   dotProduct = DOT(obj, obj2) is a Scalar object given by the Euclidean dot product of a and b. The coordinate-wise
%   multiplication and summation in this case comes from the sequence space Banach algebra.
%
%   Syntax:
%       dotProduct = DOT(obj, obj2)
%       dotProduct = obj.DOT(obj2)
%       dotProduct = DOT(obj, obj2, truncationMode) performs the multiplication using the specified truncation mode. See mtimes.m
%
%   Inputs:
%       obj - vector of length nScalar: Scalar
%       obj2 - vector of same length as obj: Scalar
%       truncationMode - specify a method of truncation: string
%
%   Outputs:
%       dotProduct - Scalar
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Aug-2018; Last revision: 08-Aug-2018

if nargin > 2
    truncationSize = varargin{1};

else
    truncationSize = 'Fixed';
end

if length(obj) ~= length(rightFactor)
    error('dot product requires vectors must be the same length')

else
    % compute pointwise products
    ptProduct(1) = mtimes(obj(1),rightFactor(1),truncationSize);
    for j = 2:length(obj)
        ptProduct(j) = mtimes(obj(j),rightFactor(j),truncationSize);
    end
    dotProduct = ptProduct.sum;
end
end %  dot

% Revision History:
%{
08-Aug-2018 - moved out of classdef file
%}
