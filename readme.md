
# SeqDE Overview
A Matlab toolbox for rigorous numerical methods based on the "sequence space" approach. Invariant sets are parameterized as coefficient sequences in some analytic function space and a proof of existence and rigorous error bounds come in the form of a-posteriori analytic error estimates.

See REFS for further details and examples.

# Classes
The basic usage revolves around several classes which efficiently perform and store important computations. Each class is briefly outlined below.

## Scalar
The scalar class gives a finite approximation representation for analytic functions of the form, $f: D^d \to \mathbb{R}^n$ where $D^d$ is the unit ball in $\mathbb{C}^d$.

### Properties
* Basis: Taylor, Fourier, Chebyshev
* Coefficient: A $d$-dimensional array of coefficients
* Dimension: The value of $d$ which is equivalent to the number of nontrivial dimensions in the coefficient array.
* NumericalClass: double or intval (requires IntLab toolbox for Matlab)
* Truncation: Integer vector of length $d$ denoting the number of coefficients in each direction.


### Hidden Properties
* Weight: Allows specification of weights for alternate $\ell_1$ weights. This should only be set to 'ones' to optimize numerical stability.


### Methods
* append
* bestfitdecay
* decay
* dot
* double
* dt
* eval
* exponent
* fixtime
* fouriertaylortimes
* imag
* intlabpoly
* intval
* intvaltimes
* inv
* leftmultiplicationoperator
* minus
* mtimes
* ndims
* norm
* plus
* real
*
* shift
* sqrt
* subdivide
* subsref
* tailratio
* uminus

### Usage

## Chart

### Properties
* Coordinate
* Truncation: {M, [N1,N2,...,Nd]}
* MaterialCoordinate{TimeSpan,[s11,s12;s21,s22;...;sd1,sd2]}
* Tau
* ErrorBound
* TimeDirection

### Hidden Properties
* Dimension
* Weight
* InitialData
* SubDivisionDepth
* SubDivisionTol
* CoefType
* BasisType %{'Taylor',etc}
* ParentId


### Methods

### Usage

## Atlas

```
function x = solve_logistic_eqn(x0, N)
x = x_0
for n = 1:N-1
    x(n+1) = dot(x, fliplr(x))/(n+1);
end
```

    
# Examples
