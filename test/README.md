To get diagnostic plots, make sure you have `Plots.jl` installed, and run

```julia
withenv(() -> include("test/runtests.jl"), "TEST_PLOTS"=>"true")
```

they will be output to `test/output`
