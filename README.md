# Algebrokeh.jl

A [AlgebraOfGraphics.jl](https://github.com/MakieOrg/AlgebraOfGraphics.jl)-style interface
to [Bokeh.jl](https://github.com/cjdoris/Bokeh.jl).

## Example

```julia
using Revise, DataFrames, Random, Bokeh, Algebrokeh

# (optional) displays the plot in the browser - omit if you are in a notebook
Bokeh.settings!(display=:browser)

# generate some example data
Random.seed!(1234)
d = DataFrame(name = repeat(["A","B","C","D","E","F"], inner=4), time=repeat([0,1,3,6], outer=6), value = rand(24))

# Create a simple scatter plot. It is automatically displayed at the REPL or in a notebook.
# From a function, call `display` on the result to show it.
plot(d, Scatter, x="@time", y="@value", color="@name")

# A more complex example which uses several different glyphs: `Scatter`, `Text` and
# `MultiLine`. The `Line` glyph does not allow `color` to be a mapping, so the `linesby()`
# function reshapes the data into the format expected by `MultiLine`.
plot(d, [Scatter, Bokeh.Text, linesby("name")], x="@time", y="@value", color="@name", text="@name")
```

![Example plot](https://raw.githubusercontent.com/cjdoris/Algebrokeh.jl/main/example.png)
