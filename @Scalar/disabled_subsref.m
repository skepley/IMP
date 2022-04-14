function varargout = subsref(obj,S)
%SUBSREF - implements subscript reference for Scalars
%
%   Subfunctions: none
%   Classes required: @Scalar
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 25-Jul-2018; Last revision: 02-Aug-2018


% disable user defined subsref until it is fixed. Bad case: type order is: '()','.','()'
[varargout{1:nargout}] = builtin('subsref', obj, S);
return

%% slice into subfields/methods via matlab builtin
if strcmp(S(1).type, '.')
    [varargout{1:nargout}] = builtin('subsref',obj,S);
    return
end

%% slice into Scalar and/or Scalar.coefficient arrays
if numel(obj) ==  1 && ~strcmp(S(1).type,'c') % slice into coefficient array of 1 Scalar
    [varargout{1:nargout}] = subsref(obj.Coefficient, S(1));

elseif strcmp(S(1).type,'()') % slice into an array of Scalars
    numDims = sum(size(obj) ~= 1); % dimensions of obj array
    if isequal(length(S(1).subs),numDims) || isequal(length(S(1).subs),1) % slice into Scalar array
        varargout{1} = builtin('subsref',obj,S(1));

    else % slice into scalar array first, then slice into coefficient array
        newSubs.type = '()';
        newSubs.subs = S(1).subs(1:numDims);
        subScalar = subsref(obj, newSubs); % pick off Scalar array from first subsref
        newSubs.subs = S(1).subs(numDims+1:end); % remaining subs
        newSubs.type = 'c'; % change call to slice coefficients
        varargout{1} = subsref(subScalar, newSubs); % pick off coefficients given by remaining subsref
    end

elseif strcmp(S(1).type,'c') % slice coefficient array for each element of Scalar array
    newSubs.type = '()';
    newSubs.subs = S(1).subs;
    if isequal(numel(obj),1)
        varargout{1} = subsref(obj, newSubs);
    else
        subOneScalar = @(scalar)subsref(scalar, newSubs);
        varargout{1} = arrayfun(subOneScalar,obj,'UniformOutput',false);
    end

else % catch all other cases
    [varargout{1:nargout}] = builtin('subsref',obj,S);
end

if length(S) > 1
    [varargout{1:nargout}] = subsref(varargout{:}, S(2:end));
end
end %  subsref

% Revision History:
%{
02-Aug-2018 - disabled subscripted reference. Some cases need to be fixed.
%}
