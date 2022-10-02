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

# You cannot use the `Line` glyph with the `color` mapping.
# For now you must explicitly group the data and use the `MultiLine` glyph - this will be
# handled automatically in the future.
gd = combine(groupby(d, :name), d->(time=[d.time], value=[d.value]))

# Create an Algebrokeh plot:
# - `data()` specifies a data source
# - `visual()` specifies a glyph to draw
# - `mapping()` specifies the mappings of glyph properties to data
# - `*` combines layers into one, to specify a combination of data, glyphs and mappings
# - `+` stacks layers on top of each other
# Then convert it to a Bokeh plot (with `draw`) and display it.
((data(d) * (visual(Scatter) + visual(Bokeh.Text)) + data(gd) * visual(MultiLine))
 * mapping(x="time", y="value", color="name", text="name")
) |> draw |> display
```

![Example plot](https://raw.githubusercontent.com/cjdoris/Algebrokeh.jl/main/example.png)
