module Algebrokeh

using Base: @kwdef
using Bokeh: Bokeh
using DataFrames: DataFrame, groupby, combine, nrow
using Tables: Tables
using DataAPI: DataAPI

export plot, linesby

include("typedefs.jl")
include("layers.jl")
include("mappings.jl")
include("draw.jl")
include("extras.jl")


# ### THEME

# Base.@kwdef mutable struct Theme
#     continuous_palette::Any = "Viridis"
#     categorical_palette::Any = "Dark2"
#     markers::Any = ["circle", "square", "triangle"]
#     legend_location::String = "right"
#     figure_opts::NamedTuple = NamedTuple()
# end

# ### MAPPING

# # ### VSTACK

# function stack(k1, k2, fields; kw...)
#     haskey(kw, k1) && error("invalid argument $k1")
#     haskey(kw, k2) && error("invalid argument $k2")
#     fields = collect(String, fields)
#     layers = Layer[]
#     x1 = Bokeh.Expr(Bokeh.Stack(fields=[]))
#     for i in 1:length(fields)
#         x2 = Bokeh.Expr(Bokeh.Stack(fields=fields[1:i]))
#         properties = Dict{Symbol,Any}()
#         for (k, v) in kw
#             if v isa AbstractVector
#                 properties[k] = v[i]
#             else
#                 properties[k] = v
#             end
#         end
#         properties[k1] = x1
#         properties[k2] = x2
#         layer = Layer(; properties)
#         push!(layers, layer)
#         x1 = x2
#     end
#     return Layers(layers)
# end

# function vstack(fields; kw...)
#     return stack(:bottom, :top, fields; kw...)
# end

# function hstack(fields; kw...)
#     return stack(:left, :right, fields; kw...)
# end

# ### DRAW

# function _get_source(data, transforms, source_cache)
#     return get!(source_cache, (data, transforms)) do 
#         if isempty(transforms)
#             if Bokeh.ismodelinstance(data, Bokeh.DataSource)
#                 return data
#             elseif Bokeh.ismodelinstance(data)
#                 error("expecting data to be a DataSource, got a $(Bokeh.modeltype(data).name)")
#             else
#                 return Bokeh.ColumnDataSource(; data)
#             end
#         else
#             source0 = _get_source(data, transforms[1:end-1], source_cache)
#             @assert Bokeh.ismodelinstance(source0)
#             Bokeh.ismodelinstance(source0, Bokeh.ColumnDataSource) || error("can only apply data transforms to ColumnDataSource, got a $(Bokeh.modeltype(source0).name)")
#             data0 = DataFrame(source0.data)
#             data = transforms[end](data0)
#             if data === data0
#                 # in-place
#                 source0.data = data
#                 return source0
#             else
#                 return Bokeh.ColumnDataSource(; data)
#             end
#         end
#     end
# end

# function histby(xcol, ncol; kw...)
#     xcol = convert(String, xcol)
#     ncol = convert(String, ncol)
#     tr = let xcol = xcol, ncol = ncol
#         datat() do df
#             return combine(groupby(df, xcol), [xcol => (x->[first(x)]) => xcol, nrow => ncol])
#         end
#     end
#     return tr * glyph(Bokeh.VBar; kw...) * mapping(xcol, ncol)
# end
# export histby

# function plotvbar(args...; x, y, dodge=nothing, kw...)
#     dodge === nothing && return plot(args..., Bokeh.VBar; x, y, kw...)
#     xcol = mapfield(x)
#     ycol = mapfield(y)
#     dodgecol = mapfield(dodge)
#     width = get(kw, :width, nothing)
#     newxcol = string(gensym("dodge#$xcol#$dodgecol"))
#     function tr(df)
#         # TODO: the list of factors should be defined somewhere extrinsic - maybe use DataFrames metadata?
#         facidxs = Dict(x=>i for (i, x) in enumerate(sort(unique(df[!, dodgecol]))))
#         nfacs = length(facidxs)
#         if nfacs > 0
#             w = something(width, 1/nfacs)
#             df[!, newxcol] = [(a, (facidxs[x] - (nfacs + 1) / 2) * w) for (a, x) in zip(df[!, xcol], df[!, dodgecol])]
#         else
#             df[!, newxcol] = []
#         end
#         return df
#     end
#     return plot(args..., Bokeh.VBar, tr; x="@$newxcol", y, kw...)
# end
# export plotvbar

end # module
