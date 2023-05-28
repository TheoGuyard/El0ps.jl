# Installation

`El0ps.jl` can be installed through Julia's [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) manager as follows:

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

The package is still in beta version and is not accessible through Julia's [general registry](https://github.com/JuliaRegistries/General) yet.
An alpha version will be released soon.
To make sure that everything went properly, you can run

```julia
pkg> test El0ps
```

The package can then be loaded into a Julia script as follows: 

```@example quickstart
using El0ps
```

!!! note 
    The package is tested against Julia `1.7` and `1.8` on `Linux` architectures.
    