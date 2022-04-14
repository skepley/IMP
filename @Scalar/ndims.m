function nDims = ndims(obj)
%DIMS - returns the number of nontrivial dimensions for the Scalar
%
%   nDims = NDIMS(obj) replaces the Matlab command ndims to return accurate dimensions for 0 or 1 dimensional Scalars.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 13-Aug-2018; Last revision: 13-Aug-2018

nDims = size(obj.Coefficient) > 1;
end % dims

% Revision History:
%{
13-Aug-2018 - moved out of classdef file
%}
