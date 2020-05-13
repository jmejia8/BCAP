# BCAP

Configure an algorithm.


## Installation


### Julia 1.0 or Later

Open the Julia REPL and press `]` to open the Pkg prompt. To add this package, use the add command:

```
pkg> add https://github.com/jmejia8/BCAP.git
```
## Quick Start


```julia
using BCAP

target_algorithm(Φ, instance, seed = 0) = begin

    # some instances solved so easy
    if instance.index % 3 == 0
        v = rand()
    else
        # Stochastic sphere function
        v = (instance.index-1)*(sum( abs.(Φ .- (0:length(Φ)-1) )  ) + 0.5rand())
    end

    # if close to the optimum put 0
    return v <= instance.optimum ? 0.0 : v
end

bounds = Array([ zeros(5) 10ones(5) ]')
parameters = Parameters(bounds, [Bool, Int, Float64, Int, Float64 ])

benchmark = [Instance(0.5i, nothing, i) for i = 1:10]

res = configure(target_algorithm, parameters, benchmark, debug = false )

```
