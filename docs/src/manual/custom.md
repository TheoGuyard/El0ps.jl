# Loss and perturbation functions

`El0ps.jl` is designed to be flexible regarding the choice of the functions $f$ and $h$.
It provides a simple way to define them.


## Functions already supported

The following functions $f$ and $h$ are already supported by the package:

| Loss / Perturbation        | Expression | Parameters
|--------------|-----|---|
| Least-Squares |  $f(\mathbf{y} ,\mathbf{A}\mathbf{x}) = \tfrac{1}{2}\|\|\mathbf{y} - \mathbf{A}\mathbf{x}\|\|_2^2$ | None |
| Logistic      |  $f(\mathbf{y} ,\mathbf{A}\mathbf{x}) = \mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y}\odot\mathbf{A}\mathbf{x}))$ | None |
| Big-M |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M)$ | $M > 0$ |
| Big-M + $\ell_1$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \alpha\|\|\mathbf{x}\|\|_1$ | $M,\alpha > 0$ |
| Big-M + $\ell_2$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \beta\|\|\mathbf{x}\|\|_2^2$ | $M,\beta > 0$ |
| $\ell_1$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1$ | $\alpha > 0$ |
| $\ell_2$-norm      |  $h(\mathbf{x}) = \beta\|\|\mathbf{x}\|\|_2^2$ | $\beta > 0$ |
| $\ell_1\ell_2$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1 + \beta\|\|\mathbf{x}\|\|_2^2$ | $\alpha,\beta > 0$ |

In the above table, the function $\mathbb{I}(\mathcal{C})$ denotes the convex indicator of the constraint $\mathcal{C}$.
The functions `f` and `h` can be instantiated as follows:
```julia
f = LeastSquares()
f = Logistic()
h = Bigm(M)
h = BigmL1norm(M, α)
h = BigmL2norm(M, β)
h = L1norm(α)
h = L2norm(β)
h = L1L2norm(α, β)
```
They all derive from the [`AbstractDatafit`](@ref) and [`AbstractPerturbation`](@ref) structs.

## Defining new loss functions

The function $f$ must verify the following hypotheses:
* $f$ is convex, proper and lower-semicontinuous
* $f$ is differentiable
* $\nabla f$ is $L$-Lipschitz

If this is fulfilled, then a new loss function can be defined as follow.
```julia
struct MyNewF <: AbstractDatafit end
```

For instance, the user can specify how it is pretty-printed with
```julia
Base.show(io::IO, f::MyNewF) = print(io, "MyNewF")
```

So that the BnB solver can run, it is also require to define the following functions:

* `lipschitz_constant(f::MyNewF, y::Vector)` : returns the value of the Lipschitz constant $L$ as a `Float64` 
* `value(f::MyNewF, y::Vector, w::Vector)` : returns the value of $f(\mathbf{y},\mathbf{w})$ as a `Float64` 
* `gradient(f::MyNewF, y::Vector, w::Vector)` : returns the value of $\nabla f(\mathbf{y},\mathbf{w})$ as a `Float64` 
* `conjugate(f::MyNewF, y::Vector, w::Vector)` : returns the value of $f^{\star}(\mathbf{y},\mathbf{w})$ as a `Float64`, where $f^{\star}(\mathbf{y},\cdot)$ denotes the convex conjugate of $f(\mathbf{y},\cdot)$
* `dual_scale!(f::MyNewF, y::Vector, w::Vector)` : in-place transforms any vector `w` into a new vector belonging to the domain of $f^{\star}(\mathbf{y},\cdot)$. No transformation is performed if the domain of $f^{\star}(\mathbf{y},\cdot)$ is the whole space.

Then, the BnB is able to handle losses that are instances of `MyNewF` on his own.

## Defining new perturbation terms

The function $h$ must verify the following hypotheses:
* $h$ splits, i.e., $h(\mathbf{x}) = \sum_{i} h_i(x_i)$
* The splitting terms are proper, convex and lower-semicontinuous
* The splitting terms are equals, i.e., $h_i=h_j$ for all $i \neq j$
* The splitting terms are even
* The splitting terms are coercive
* The splitting terms verify $h_i(x_i) \geq h_i(0) = 0$, for all $i$

The above hypotheses are for instance verifies for norms, bound constraints, ...
If they are fulfilled, then a new perturbation term can be defined as follow.
```julia
struct MyNewH <: AbstractPerturbation end
```

For instance, the user can specify how it is pretty-printed with
```julia
Base.show(io::IO, H::MyNewH) = print(io, "MyNewH")
```

So that the BnB solver can run, it is also require to define the following functions:
* `value_1d(h::MyNewH, x::Float64)`: returns the value of $h_i(x_i)$ as a `Float64` 
* `conjugate_1d(h::MyNewH, x::Float64)`: returns the value of the confjugate function $h_i^{\star}(x_i)$ as a `Float64` 
* `prox_1d(h::MyNewH, x::Float64, η::Float64)`: returns the proximity operator $\mathrm{prox}_{\eta h_i(\cdot)}(x)$ as a `Float64`
* `dual_scale!(h::MyNewH, A::Matrix, u::Vector, λ::Float64)`: in-place transforms any vector `u` into a new vector belonging to the domain of $f^{\star}(\mathbf{y},\tfrac{1}{\lambda}\mathbf{A}^{\top}\cdot)$. Returns the value of $\mathbf{A}^{\top}\mathbf{u}$

Then, the BnB is able to handle losses that are instances of `MyNewH` on his own.

## Examples

### Loss function

The following portion of code shows how to implement a quadratic loss function $f(\mathbf{y},\mathbf{w}) = \tfrac{1}{2}\mathbf{w}^{\top}\mathbf{Q}\mathbf{w} + \mathbf{w}^{\top}\mathbf{y}$.
In this case, one has
* $L = \sigma_{\max} (\mathbf{Q})$, where $\sigma_{\max}(\cdot)$ is the largest singular of some matrix.
* $\nabla f(\mathbf{y},\mathbf{w}) = \mathbf{Qw + y}$
* $f^{\star}(\mathbf{y},\mathbf{w}) = \mathbf{Q}^{-1}\mathbf{(w - y)}$

and the domain of $f^{\star}$ is the whole space.

```julia
# Definition of the new struct 
struct QuadraticLoss <: AbstractDatafit 
    Q::Matrix
    QuadraticLoss(Q::Matrix) = new(Q)
end

# Definition of the loss operators
lipschitz_constant(f::QuadraticLoss, y::Vector) = maximum(svdvals(f.Q))
value(f::QuadraticLoss, y::Vector, w::Vector) = 0.5 * (w' * f.Q * w) + w' * y
gradient(f::QuadraticLoss, y::Vector, w::Vector) = f.Q * w + y
conjugate(f::QuadraticLoss, y::Vector, w::Vector) = inv(Q) * (w - y)
dual_scale!(f::QuadraticLoss, y::Vector, w::Vector) = nothing
```

Note that computations can be easily casev by computing and storing the value of $L$ and $\mathbf{Q}^{-1}$ once for all when initializing the `QuadraticLoss` struct.

### Perturbation term

The following portion of code shows how to implement a p-norm perturbation term $h(\mathbf{x}) = \sum_{i}h_i(x_i)$
 with $h_i(x_i) = \tfrac{1}{p}|x_i|^p$ for all $i$.
In this case, one has
* $h_i^{\star}(u) = \tfrac{1}{q}|u|^q$
* $\mathrm{prox}_{\eta h_i(\cdot)}(x) = \mathrm{sign}(x) \times \mathrm{root}(\eta r^{p-1}+r-|x|)$

where $q$ is such that $\tfrac{1}{p} + \tfrac{1}{q} = 1$ and where $\mathrm{root}(P(r))$ returns a root of the polynomial $P(r)$.

```julia
# For the prox operation
using Roots

# Definition of the new struct 
struct LpNorm <: AbstractPerturbation
    p::Int
    q::Int
    LpNorm(p::Int) = new(p, p/(p-1))
end

# Definition of the loss operators
value_1d(h::LpNorm, x::Float64) = abs(x)^(h.p) / h.p
conjugate_1d(h::LpNorm, x::Float64) = abs(x)^(h.q) / h.q
prox_1d(h::LpNorm, x::Float64, η::Float64) = sign(x) * find_zero(r -> ηr^(h.p-1) + r - abs(x), 0) 
dual_scale!(h::MyNewH, A::Matrix, u::Vector, λ::Float64) = A' * u
```

Above, the `Roots` package is used in the `prox_1d` to solve the polynomial equation.