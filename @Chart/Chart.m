classdef Chart < handle
    %CHART - A class for representing a complex analytic function, f: C^d ---> C^n.
    %
    %   CHART() - A more detailed description of the class
    %
    %   CHART constructor syntax:
    %       CHARTObj = $NAME()
    %
    %   CHART properties:
    %       Coordinate - A length n vector of Scalars whose domain is a subset of [-1,1]^d.
    %       Dimension - [d,n] with d <= n are the domain and codomain dimension for the Chart.
    %       Truncation - [N1,N2,...,Nd] is a vector identifying the truncation space in which to consider the coordinates of this Chart.
    %       TimeSpan - [t0, tf] global coordinates for 1st dimension variable which is often time-like. The "globalness" here is due to the correspondence 
    %           between [t0, tf] and the interval [0,1] which are the local coordinates for this chart. If tf < t0 (i.e. integration is backwards in time), 
    %           the local parameterization is still on [0,1] with P(0,:) always pointing to the initial data. 
    %       Tau - nonzero real number which is negative for backward time integration and positive for forward time integration.
    %
    %   CHART methods:
    %       Method 1 - description
    %       Method 2 - description
    %
    %   Examples:
    %       Line 1 of example
    %       Line 2 of example
    %       Line 3 of example
    %
    %   Subclasses: none
    %   Superclasses: none
    %   Other classes required: Scalar
    %   Other m-files required: none
    %   MAT-files required: none
    %
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 21-Aug-2018; Last revision: 19-Apr-2019
    %
    %   ToDo:
    %   item1 -
    %   item2 -
    
    
    properties
        Coordinate; % length n vector of Scalars 
        Basis; % A cell array of analytic base functions for each dimension
        Dimension; % [d, n] 
        Truncation; % [N1,...,Nd]
        TimeSpan; % [t0, tf] 
        SpatialSpan; % global coordinates with respect to [-1,1]^d for dimensions 2,3,...,d which are often space-like.
        Tau; % |t1 - t0|
        Crash = false; % Generic property to flag charts which which aren't well conditioned for some context specific reason.
        ErrorBound = []; % rigorous error bound for validated charts or empty if not validated.
        NumericalClass; % double or intval
    end
    
    properties(Hidden = 1)
        Weight = 'ones'; % ell^1 space weights
        InitialData; % length n vector of Scalars of dimension d
        SubDivisionDepth; % recursion depth for use by subdivision algorithms
        InitialError = 0; % initial ell_1 error
        StepError = 0; % Validation error for this timestep (ErrorBound - InitialError)
        ParentHandle; % handle pointing to the parent chart if it exists
        % PtWiseNorm; % scalars of the form (|a_j| + r_j)
        IsValid;
        ClassConstructor; % A function handle to the constructor to create more objects of this same subclass
        Generation = 0;
    end
    
end %end class

% Revision History:
%{
7-Mar-2019 - Re-factored for use with the new Atlas class. Large scale overhaul of this entire class to include a lot of bug fixes and code 
    class/method refactoring. 
22-Mar-2019 - Removed Flow Direction property and refactored class to allow negative values for tau. 
19-Apr-2019 - Changed initial Error Bound to be an empty array instead of 0. 
%}

