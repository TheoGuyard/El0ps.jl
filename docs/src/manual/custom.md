# Loss and perturbation functions

`El0ps.jl` is designed to be flexible regarding the choice of the functions $f$ and $h$.

## Functions already supported

The following functions $f$ and $h$ are already supported by the package:

| Loss / Perturbation        | Expression | Parameters
|--------------|-----|---|
| Least-Squares |  $f(\mathbf{A}\mathbf{x}) = \tfrac{1}{2}\|\|\mathbf{y} - \mathbf{A}\mathbf{x}\|\|_2^2$ | None |
| Logistic      |  $f(\mathbf{A}\mathbf{x}) = \mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y}\odot\mathbf{A}\mathbf{x}))$ | None |
| Big-M |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M)$ | $M > 0$ |
| Big-M + $\ell_1$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \alpha\|\|\mathbf{x}\|\|_1$ | $M,\alpha > 0$ |
| Big-M + $\ell_2$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \beta\|\|\mathbf{x}\|\|_2^2$ | $M,\beta > 0$ |
| $\ell_1$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1$ | $\alpha > 0$ |
| $\ell_2$-norm      |  $h(\mathbf{x}) = \beta\|\|\mathbf{x}\|\|_2^2$ | $\beta > 0$ |
| $\ell_1\ell_2$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1 + \beta\|\|\mathbf{x}\|\|_2^2$ | $\alpha,\beta > 0$ |

In the above table, the function $\mathbb{I}(\mathcal{C})$ denotes the convex indicator of the constraint $\mathcal{C}$.
They can simply be instantiated as follows:
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
They all derive from the [`AbstractDatafit`](@ref) and [`AbstractPenalty`](@ref) structs.

## Defining new functions

In addition, we provide a simple way to define new functions $f$ and $g$.

### Loss function

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
Base.show(io::IO, F::LeastSquares) = print(io, "MyNewF")
```

So that the BnB solver can run, it is also require to define the following functions:

* `lipschitz_constant(f::MyNewF, y::Vector)` : returns the value of the Lipschitz constant $L$ as a `Float64` 
* `value(f::MyNewF, y::Vector, w::Vector)` : returns the value of $f(\mathbf{y},\mathbf{w})$ as a `Float64` 
* `gradient(f::MyNewF, y::Vector, w::Vector)` : returns the value of $\nabla f(\mathbf{y},\mathbf{w})$ as a `Float64` 
* `conjugate(f::MyNewF, y::Vector, w::Vector)` : returns the value of $f^{\star}(\mathbf{y},\mathbf{w})$ as a `Float64`, where $f^{\star}(\mathbf{y},\cdot)$ denotes the convex conjugate of $f(\mathbf{y},\cdot)$
* `dual_scale!(f::MyNewF, y::Vector, w::Vector)` : in-place transforms any vector `w` into a new vector belonging to the domain of $f^{\star}(\mathbf{y},\cdot)$. No transformation is performed if the domain of $f^{\star}(\mathbf{y},\cdot)$ is the whole space.

Then, the BnB is able to handle losses that are instances of `MyNewF` on his own.

### Perturbation term

TODO


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
lipschitz_constant(F::QuadraticLoss, y::Vector) = maximum(svdvals(F.Q))
value(F::QuadraticLoss, y::Vector, w::Vector) = 0.5 * (w' * F.Q * w) + w' * y
gradient(F::QuadraticLoss, y::Vector, w::Vector) = F.Q * w + y
conjugate(F::QuadraticLoss, y::Vector, w::Vector) = inv(Q) * (w - y)
dual_scale!(F::QuadraticLoss, y::Vector, w::Vector) = nothing
```

Note that computations can be easily casev by computing and storing the value of $L$ and $\mathbf{Q}^{-1}$ once for all when initializing the `QuadraticLoss` struct.