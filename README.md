# ORNL_2024_Lectures

## Materials from 2024 Lecture Series at ORNL, "Ideas behind the Sunny.jl package"

This repository contains the lecture slides (in the `data` directory) as well as the code examples used during the lectures (in the `scripts` directory). These examples
were rather informal. For more polished illustrations of how to use [Sunny](https://github.com/SunnySuite/Sunny.jl), please see the [official docs](https://sunnysuite.github.io/Sunny.jl/dev/).

For guidance on setting up a Julia environment and installing Sunny,
please see the documentation [here](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia). Scripts in this repository that have names beginning with either `01` or `02` should work in a regular Julia environment.
In these cases, it will not be necessary to execute the following two lines, which are found at the top of each example.
```
using DrWatson
@quickactivate "ORNL_2024_LECTURES"
```
The "entangled units" example requires installing a development branch of Sunny. This can be achieved
most easily by installing the DrWatson package and executing the lines above. Instructions
for installing DrWatson are included below.


## Instructions for use
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ORNL_2024_Lectures

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ORNL_2024_LECTURES"
```
which auto-activate the project and enable local path handling from DrWatson.
