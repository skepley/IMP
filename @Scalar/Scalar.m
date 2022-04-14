classdef Scalar
    %SCALAR - SeqDE class for representing analytic scalar functions.
    %
    %   The Scalar class is a finite approximation for an analytic scalar of the form, f: D^d --> R^n where D^d is
    %   a d-dimensional polydisc in C^d. Scalar representations are specified as double or intval coefficients
    %   with respect to a Taylor, Fourier, or Chebyshev series expansion.
    %
    %   SCALAR constructor syntax:
    %       ScalarObj = SCALAR(C, basis) constructs a Scalar whose coefficient sequence is specified by C with respect to the basis
    %
    %       ScalarObj = SCALAR(C, basis, [N1,...,Nd]) embed the Scalar in the d-dimensional truncation space by padding with zeros as necessary.
    %
    %   SCALAR properties:
    %       coefficient - array of coefficients
    %       Dimension - dimension of the polydisc
    %       Truncation - vector of truncation degrees for each dimension
    %       Weight - vector of weights for sequence space (default is [1,...,1])
    %       NumericalClass - string specifying double, interval, or Scalar coefficients
    %
    %   SCALAR methods:
    %       intval - convert to interval coefficients
    %
    %   Examples:
    %       A = SCALAR([1,2,3], 'Taylor') returns a Scalar representation of the univariate polynomial f(s) = 1 + 2*s + 3*s^2.
    %
    %       A = SCALAR(magic(3), 'Taylor') returns a Scalar representation of the multivariate polynomial
    %           f(s,t) = (8 + s + 6s^2) + (3 + 5s + 7s^2)t + (4 + 9s + 2s^2)t^2
    %
    %       A = SCALAR([1,2,3], 'Basis', [3,3,3]) returns the Scalar representation for the function, f: D ---> R, embedded
    %           as a function on the domain, D^3. This corresponds to padding the polynomial with zeros in the additional variables i.e.
    %           f(s1,s2,s3) = 1 + 2s + 3s^2
    %
    %   Subclasses: none
    %   Superclasses: none
    %   Other classes required: none
    %   Other m-files required: INTLAB toolbox
    %   MAT-files required: none
    %
    %   See also: @Chart, @Atlas
    
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 08-Jun-2016; Last revision: 18-Aug-2018
    %
    %   ToDo:
    %   Merge and simplify evaluation, fixtime, derivative, and truncation methods/function calls.
    %   Rigorous evaluation of sin/cos/exp, inverses, and fractional powers.
    %   Finish vectorization of all methods.
    %   Finish high precision inteval FFT and reimplement FFT-based fast convolution.
    %   Add polynomial evaluation over the Banach algebra.
    %   Fix shift method to use new subscript reference.
    
    
    %% -------------------- Properties --------------------
    properties
        Basis; % Taylor, Fourier, Chebyshev: string
        Coefficient; % N1-by-N2-by-...-Nd array: double or intval.
        Dimension; % d: integer
        NumericalClass; % double, intval, Scalar: string
        Truncation; % [N1,...,Nd]: integer
    end
    
    properties(Hidden = 1)
        %         Weight = 'ones'; % [nu_1,...,nu_d]: positive double
    end
    
    %% -------------------- Methods --------------------
    methods
        function obj = Scalar(coefData, basis, varargin)
            %SCALAR - class constructor
            
            if nargin > 0
                narginchk(2, 3)
                
                if isa(coefData, 'cell') % vectorized constructor
                    arrayDim = size(coefData);
                    coefData = reshape(coefData, [], 1); % flatten data to a column vector
                    nCell = length(coefData);
                    obj(nCell) = Scalar(coefData{nCell}, basis, varargin{:}); % preallocation and assigning last index
                    
                    % loop over remaining cells
                    for iCell = 1:nCell-1
                        obj(iCell) = Scalar(coefData{iCell}, basis, varargin{:});
                    end
                    obj = reshape(obj, arrayDim); % reshape to original size
                    return
                    
                elseif isa(coefData, 'Scalar') % coefficient specifed as Scalar
                    if length(coefData) > 1 % coefData is a vector of Scalars
                        obj = coefData;
                        for iScalar = 1:numel(obj)
                            obj(iScalar) = Scalar(obj(iScalar), basis, varargin{:});
                        end
                        return
                    else
                        if ~strcmp(basis, coefData.Basis)
                            error('Basis conversion not implemented')
                            
                        elseif nargin == 2
                            obj = coefData; % return a deep copy
                            
                        else % new embedding specified
                            obj = coefData; % copy of the Scalar
                            obj.Truncation = varargin{1};
                            obj.Dimension = length(obj.Truncation);
                            obj.Coefficient = Scalar.embed(coefData.Coefficient, obj.Truncation);
                        end
                        return
                    end
                end
                
                % set Scalar properties
                switch class(coefData)
                    case 'polynom'   % coefficient is an IntLab polynomial
                        if ~strcmp(basis, 'Taylor')
                            error('Intlab polynom class is only compatiable with Taylor basis')
                        end
                        coefTruncation = 1 + max(coefData.e);
                        flippedCoefArray = squeeze(reshape(coefData.c, coefTruncation));
                        coefArray = Scalar.flipcoef(flippedCoefArray, 1:ndims(flippedCoefArray));
                        obj.NumericalClass = class(coefArray);
                        
                        if nargin == 3
                            obj.Truncation = varargin{1};
                            obj.Dimension = length(obj.Truncation);
                            obj.Coefficient = Scalar.embed(coefArray, obj.Truncation); % pad coefficients to correct size
                        else
                            dims = size(coefArray);
                            trueDims = dims(dims > 1); % remove singleton dimensions
                            obj.Dimension = length(trueDims);
                            if obj.Dimension > 0
                                obj.Truncation = trueDims;
                            else
                                obj.Truncation = 1; % truncation = 1 for 0-dimensional (constant) Scalar
                            end
                            obj.Coefficient = coefArray;
                        end
                        obj.Basis = basis;
                        
                        
                    otherwise % coefficient specifed as double or intval array
                        obj.NumericalClass = class(coefData); % intval or double
                        if nargin == 3
                            obj.Truncation = varargin{1}; % specify truncation as vararg.
                            obj.Dimension = length(obj.Truncation);
                            obj.Coefficient = Scalar.embed(coefData, obj.Truncation); % pad coefficients to correct size
                        else
                            coefData = squeeze(coefData);
                            dims = size(coefData);
                            trueDims = dims(dims > 1); % remove singleton dimensions
                            obj.Dimension = length(trueDims);
                            if obj.Dimension > 0
                                obj.Truncation = trueDims;
                            else
                                obj.Truncation = 1; % truncation = 1 for 0-dimensional (constant) Scalar
                            end
                            obj.Coefficient = coefData;
                        end
                        
                        % write basis
                        if isa(basis, 'cell')
                            obj.Basis = basis;
                        elseif isequal(obj.Dimension, 0)
                            obj.Basis = {basis}; % handle 0 dimensional initial data
                        else
                            obj.Basis = cell(1,obj.Dimension);
                            for j = 1:obj.Dimension
                                obj.Basis{j} = basis;
                            end
                        end
                end % switch class(coefData)
            end % if nargin > 0
        end %  class constructor
        
        
        %% -------------------- METHOD REPAIR SHOP --------------------
        % The following methods need to be updated and thoroughly checked before using.
        
        function int_ds = intds(obj, varargin)
            % evaluate definite or indefinite integral with respect to spatial variable
            if obj.Dimension == 2
                C = repmat(1./(1:obj.Truncation(2)),obj.Truncation(1),1);
                int_ds = Scalar([zeros(obj.Truncation(1),1), basis, C.*obj.Coefficient]); % indefinite integral
                if (nargin > 1) % compute integral obj ds on [a,b]
                    bounds = varargin{1};
                    int_ds = int_ds.fix_space(bounds(2)) - int_ds.fix_space(bounds(1)); % evaluation definite integral
                end
            else
                error('ds not implemented for this dimension')
            end
        end
        
        
        % THIS NEEDS TO BE IMPLEMENTED USING INTLAB POLYEVAL AND MERGED INTO THE EVAL METHOD.
        function evalObj = intvalEval(obj,s,t)
            % s is a vector of intervals, t is a single scalar.
            X = obj.fixtime(t);
            F = polynom(flip(X.Coefficient),'s');
            evalObj = polyval(F,s);
        end
        
        % THIS NEEDS TO BE MERGED INTO THE EVAL METHOD.
        function evalObj = fixSpace(obj,s)
            % collapse Taylor series onto fixed space evaluation, s
            space_vector = bsxfun(@power,s,(0:size(obj.Coefficient,2)-1)');
            evalObj = obj.Coefficient*space_vector;
        end
        
        
        
        
        
        %% -------------------- METHOD GRAVEYARD --------------------
        % The following methods should not be called but are retained for a while to make sure they don't break any old code.
        
        function newObj = homog(obj)
            % convert Scalar coefs to double or intval coefs
            warning('homog is deprecated and Scalar NumercalClass is not currently allowed')
        end
        
        
        
        function truncate(obj,modes)
            % truncate Scalar to specified number of modes
            
            warning('truncate method is deprecated. Use the static embed method or the Scalar class constructor')
            if ~isequal(obj.Dimension,2)
                error('not implemented for dimension other than 2 yet')
            end
            
            if isequal(length(obj),1) % a single Scalar
                newTruncation = [min([obj.Truncation(1),modes(1)]),min([obj.Truncation(2),modes(2)])];
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Coefficient = obj.Coefficient(1:newTruncation(1),1:newTruncation(2));
                obj.Truncation = newTruncation;
            else
                for j = 1:length(obj)
                    truncate(obj(j),modes)
                end
            end
        end
        
        
        function columnObj = col(obj)
            % returns coefficients of BAsclar as a column vector (double) under the canonical isomorphism
            
            warning('col is deprecated and should be replaced by the exponent method')
            columnObj = reshape(obj.Coefficient,[],1);
        end
        
        
        function dObj_ds = ds(obj)
            % compute spatial derivative
            
            warning('ds is deprecated. Use dt method instead.')
            switch obj.Dimension
                case 1
                    C = 1:obj.Truncation-1;
                    dObj_ds = Scalar(C.*obj.Coefficient(2:end), obj.Basis);
                case 2
                    C = repmat((1:obj.Truncation(2)-1),obj.Truncation(1),1);
                    dObj_ds = Scalar(C.*obj.Coefficient(:,2:obj.Truncation(2)), obj.Basis);
                otherwise
                    error('ds not implemented for this dimension')
            end
        end
        
        
        
    end %  methods
    
    %% STATIC METHODS
    methods(Static)
        function embedCoef = embed(coefData, truncation)
            % pads coefficient with zeros to achieve specified truncation size
            
            if isequal(truncation, 1) % constant Scalar
                embedCoef = coefData;
                return
            end
            
            % coefficient subscript reference
            truncDim = length(truncation); % true dimension of the coefficient array
            coefDim = length(size(coefData));
            dimDeficit = truncDim - coefDim; % number of dimensions missing from the the coefficient data
            
            if dimDeficit >= 0
                % Handle the case that the truncation vector is longer than the number of dimensions in coefData. This embeds the coefficient data
                % into a space whose dimension is specified by the truncation vector by padding additional dimensions with initial size = 1.
                
                if isa(coefData, 'double') % initialize coefficient array of the correct size
                    embedCoef = zeros(truncation);
                elseif isa(coefData, 'intval')
                    embedCoef = intval(zeros(truncation));
                end
                coefTrunc = min(size(coefData), truncation(1 + truncDim - coefDim : end));
                coefSubs = arrayfun(@(dim)1:coefTrunc(dim), 1:length(coefTrunc), 'UniformOutput', false);
                
                % embedding subscript reference
                missingDim = ones(1, length(truncation) - coefDim); % padded dimensions
                subTruncation = [missingDim, coefTrunc]; % padded dimensions append indices to the left
                padSubs = arrayfun(@(dim)1:subTruncation(dim), 1:length(subTruncation), 'UniformOutput', false);
                
                % set coefficients
                embedCoef(padSubs{:}) = coefData(coefSubs{:});
            elseif isequal(truncDim, 1) && isequal(dimDeficit, -1)
                % deficit is negative should happen only if the dimension is 1 so the coefficient data is a row or column vector which MATLAB
                % thinks has dimension 2.
                
                if isa(coefData, 'double') % initialize coefficient array of the correct size
                    embedCoef = zeros(1, truncation);
                elseif isa(coefData, 'intval')
                    embedCoef = intval(zeros(1, truncation));
                end
                embedCoef(1:min(truncation, length(coefData))) = coefData(1:min(truncation, length(coefData))); % slice coefficient data into array
                
                if iscolumn(coefData)
                    embedCoef = embedCoef'; % reshape to match coefficient data
                end
            else
                error('Number of truncation dimensions is underspecified')
            end
        end % embed
        
        
        function randiScalar = randi(varargin)
            % quickly generate a Taylor Scalar for testing Scalar functions and methods
            randiScalar = Scalar(randi(varargin{:}), 'Taylor');
        end % randi
        
        
        function zeroScalar = zeros(basis, varargin)
            % initialize a zero Scalar
            zeroScalar = Scalar(zeros(varargin{:}), basis);
        end
        
        function newCoef = flipcoef(coefData, idx)
            % flip a coefficient array along each specified index
            if isequal(length(idx),1)
                newCoef = flip(coefData, idx(1)); % flip dimension
            else
                newCoef = Scalar.flipcoef(flip(coefData, idx(end)), idx(1:end-1)); % flip and recurse
            end
        end
        
        function truncSubs = trunc2subs(varargin)
            % given a truncation row vector [N1,...,Nd] returns a cell array of indices {1:N1,1:N2,...,1:Nd}
            % given vectors [n1,...,nd], [N1,...,Nd] returns a cell array {n1:N1,...,nd:Nd}
            
            if nargin == 1
                truncHigh = varargin{1};
                truncLow = ones(size(truncHigh));
            else
                [truncLow, truncHigh] = varargin{:};
            end
            truncSubs = arrayfun(@(L,H)L:H,truncLow, truncHigh, 'UniformOutput', false);
        end
        
        function dimension = dims(coefArray)
            % return the true dimension of coefficient array
            fullDim = size(coefArray);
            trueDim = fullDim(fullDim > 1); % remove singleton dimensions
            if isempty(trueDim)
                dimension = 0;
            else
                dimension = length(trueDim);
            end
        end
        
        function [a,b] = commonsize(a,b)
            % resizes two arrays to their smallest common embedding by padding with zeros.
            szA = Scalar.trunc2subs(size(a)+1, size(b));
            szB = Scalar.trunc2subs(size(b)+1, size(a));
            a(szA{:}) = 0;
            b(szB{:}) = 0;
        end
        
    end %  static methods
end %  classdef

% Revision History:
%{
11-Jul-2017 - Support for interval and Scalar coefficients added.
15-Aug-2017 - Reverted FFT based convolution to a classical algorithm. FFT is faster but numerically unstable, especially for intval
    coefficients.
08-Aug-2018 - Class completely overhauled and renamed from BAScalar to Scalar. Numerous improvements including:
   full class code refactorization and organization of class folder
   class inheritance changed from handle to value
   support for additional bases (Fourier or Chebyshev)
   Intlab polynom dependency removed
   support for higher dimensions
   streamlined truncation methods
   subscript reference methods
   numerous efficiency gains, bug fixes, and improvements in data type consistency
20-Aug-2018 - Completely rewritten class constructor. Important changes include:
    suppose for arbitrary dimensional Scalars
    1-dimensional Scalar coefficients are strictly column vectors
    bug fixes for embedding Scalars into truncation spaces
    addition of overloaded zeros, randi methods
16-Jan-2019 - Fixed a bug in the embed static method when initializing 1-dimensional Scalars and specifying a truncation size. Bug fixes for
    decay, bestfitdecay, and embed methods. Details in their respective files.
20-Feb-2020 - Fixed a bug where 0 dimensional Scalars would not have a basis defined.
%}


