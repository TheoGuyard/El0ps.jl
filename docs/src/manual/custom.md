# Loss and penalty functions

This package is designed to be flexible regarding the choice of the functions $f$ and $h$.
It provides a simple way to define them.


## Default support

The package already implements by default the following functions:

| Loss / Penalty        | Expression | Parameters
|:--------------|:-----|:---|
| Least-Squares |  $f(\mathbf{A}\mathbf{x}) = \tfrac{1}{2}\|\|\mathbf{y} - \mathbf{A}\mathbf{x}\|\|_2^2$ | Vector $\mathbf{y}$ |
| Logistic      |  $f(\mathbf{A}\mathbf{x}) = \mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y}\odot\mathbf{A}\mathbf{x}))$ | Vector $\mathbf{y}$ |
| Big-M |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M)$ | Scalar $M > 0$ |
| Big-M + $\ell_1$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \alpha\|\|\mathbf{x}\|\|_1$ | Scalars $M,\alpha > 0$ |
| Big-M + $\ell_2$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \beta\|\|\mathbf{x}\|\|_2^2$ | Scalars $M,\beta > 0$ |
| $\ell_1$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1$ | Scalar $\alpha > 0$ |
| $\ell_2$-norm      |  $h(\mathbf{x}) = \beta\|\|\mathbf{x}\|\|_2^2$ | Scalar $\beta > 0$ |
| $\ell_1\ell_2$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1 + \beta\|\|\mathbf{x}\|\|_2^2$ | Scalars $\alpha,\beta > 0$ |

In the above table, $\mathbb{I}(\mathcal{C})$ denotes the convex indicator of the constraint $\mathcal{C}$ and $\odot$ denotes the Hadamard product.
They can be instantiated as follows:

```julia
f = LeastSquares(y)
f = Logistic(y)
h = Bigm(M)
h = BigmL1norm(M, α)
h = BigmL2norm(M, β)
h = L1norm(α)
h = L2norm(β)
h = L1L2norm(α, β)
```

The function `f` and `h` respectively derive from the [`AbstractDatafit`](@ref) and [`AbstractPenalty`](@ref) structures.

## Defining new loss functions

To define new functions $f$, one must verify the following hypotheses:
* the function is convex, proper and lower-semicontinuous
* the function is differentiable
* the function gradient is Lipschitz-continuous

If so, then a new function $f$ can be defined as follow.

```julia
struct MyNewF <: AbstractDatafit end
```

For instance, the user can specify how it is pretty-printed with

```julia
Base.show(io::IO, f::MyNewF) = print(io, "MyNewF")
```

So that the [`BnbSolver`](@ref) can run, it is also require to define the following functions:
* `dim_input(f::LeastSquares)` : the dimension of the function input as an `Int`
* `lipschitz_constant(f::MyNewF)` : returns the value of the Lipschitz constant of the gradient of `f` as a `Float64` 
* `value(f::MyNewF, w::Vector)` : returns the value of `f` at `w` as a `Float64` 
* `gradient(f::MyNewF, w::Vector)` : returns the gradient of `f` at `w` as a `Vector` 
* `conjugate(f::MyNewF, w::Vector)` : returns the value of the convex conjugate of `f` at `w` as a `Float64` 

Once these functions have been overloaded, the [`BnbSolver`](@ref) is able to handle losses that are instances of `MyNewF` on his own.

## Defining new penalty terms

To define new functions $h$, one must verify the following hypotheses:
* the function splits
* the splitting terms are proper, convex and lower-semicontinuous
* the splitting terms are equals
* the splitting terms are even
* the splitting terms are coercive
* the splitting terms are non-negative and equal zero at zero

The above hypotheses are verified for norms or for bound constraints, among many others.
If they are fulfilled, a new penalty term can be defined as follow.

```julia
struct MyNewH <: AbstractPenalty end
```

For instance, the user can specify how it is pretty-printed with
```julia
Base.show(io::IO, H::MyNewH) = print(io, "MyNewH")
```

So that the [`BnbSolver`](@ref) can run, it is also require to define the following functions:
* `value_1d(h::MyNewH, x::Float64)`: returns the value of a splitting term evaluated at `x` as a `Float64` 
* `conjugate_1d(h::MyNewH, x::Float64)`: returns the value of the conjugate function of a splitting term evaluated at `x` as a `Float64` 
* `prox_1d(h::MyNewH, x::Float64, η::Float64)`: returns the proximity operator of `η` times a splitting term evaluated at `x` as a `Float64` 

Once these functions have been overloaded, the [`BnbSolver`](@ref) is able to handle penalty terms that are instances of `MyNewH` on his own.

## Examples

### Loss function

The following portion of code shows how to implement a quadratic loss function $f(\mathbf{w}) = \tfrac{1}{2}\mathbf{w}^{\top}\mathbf{Q}\mathbf{w} + \mathbf{w}^{\top}\mathbf{y}$.

```julia
# Definition of the new struct 
struct QuadLoss <: AbstractDatafit 
    Q::Matrix
    y::Vector
    QuadLoss(Q::Matrix, y::Vector) = new(Q, y)
end

# Definition of the operators
dim_input(f::QuadLoss) = length(f.y)
lipschitz_constant(f::QuadLoss) = maximum(svdvals(f.Q))
value(f::QuadLoss, w::Vector) = 0.5 * (w' * f.Q * w) + w' * f.y
gradient(f::QuadLoss, w::Vector) = f.Q * w + f.y
conjugate(f::QuadLoss, w::Vector) = inv(Q) * (w - f.y)
```

Note that computations can be saved by computing and storing the value of the Lipschitz constant and of $\mathbf{Q}^{-1}$ once for all when initializing the `QuadLoss` structure.

### Penalty term

The following portion of code shows how to implement a p-norm penalty term $h(\mathbf{x}) = \tfrac{1}{p}\|\mathbf{x}\|^p_p$ for some integer $p \geq 1$.

```julia
# Required for the prox operation
using Roots

# Definition of the new struct 
struct LpNorm <: AbstractPenalty
    p::Int
    q::Int
    LpNorm(p::Int) = new(p, p/(p-1))
end

# Definition of the operators
value_1d(h::LpNorm, x::Float64) = abs(x)^(h.p) / h.p
conjugate_1d(h::LpNorm, x::Float64) = abs(x)^(h.q) / h.q
prox_1d(h::LpNorm, x::Float64, η::Float64) = sign(x) * find_zero(r -> ηr^(h.p-1) + r - abs(x), 0) 
```

See the [prox-website](http://proximity-operator.net) for more details.