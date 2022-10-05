module Algebrokeh

using Base: @kwdef
using Bokeh: Bokeh
using DataFrames: DataFrame, groupby, combine, nrow
using Tables: istable

export plot


### TYPES

@enum DataType NUMBER_DATA FACTOR_DATA

@enum MappingType COLOR_MAP MARKER_MAP HATCH_PATTERN_MAP DATA_MAP

"""
Holds information about a collection of data, such as a column in a DataSource.
"""
@kwdef struct DataInfo
    datatype::Union{Nothing,DataType} = nothing
    factors::Union{Nothing,Vector} = nothing
end

"""
A DataSource, plus information about the data in its columns.
"""
struct Data
    source::Union{Bokeh.ModelInstance}
    columns::Dict{String,DataInfo}
end

"""
A specifier for a field. Can also be a collection of 2 or 3 fields, to treat them as a
hierarchical factor.
"""
struct Field
    name::Union{String,Tuple{String,String},Tuple{String,String,String}}
end

"""
A mapping value, which defines some connection between a field of data and a visual
property, such as color, marker or location.
"""
@kwdef struct Mapping
    name::Union{Nothing,String}
    type::Union{Nothing,MappingType}
    field::Field
    label::Any = nothing
    transforms::Vector{Bokeh.ModelInstance} = Bokeh.ModelInstance[]
    datainfo::DataInfo = DataInfo()
end

"""
A plotting layer, consisting of some combination of data, transforms to apply to the data,
the glyph to use to plot the data, and properties of that glyph including mappings from
data columns to visual attributes.

Any properties which are not [`Mapping`](@ref) are passed through to Bokeh unchanged.
"""
@kwdef struct Layer
    data::Union{Nothing,Data} = nothing
    transforms::Vector{Any} = []
    glyph::Union{Nothing,Bokeh.ModelType} = nothing
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

"""
A collection of layers, stacked one on top of the next. This is returned by [`plot`](@ref).
"""
struct Layers
    layers::Vector{Layer}
end

const Factor = Union{AbstractString,NTuple{2,AbstractString},NTuple{3,AbstractString}}


### LAYER

function Base.:(*)(x::Layer, y::Layer)
    data = y.data !== nothing ? y.data : x.data
    transforms = y.data !== nothing ? y.transforms : vcat(x.transforms, y.transforms)
    glyph = y.glyph !== nothing ? y.glyph : x.glyph
    properties = merge(x.properties, y.properties)
    return Layer(; data, transforms, glyph, properties)
end


### LAYERS

Base.zero(::Type{Layers}) = Layers(Layer[])

Base.one(::Type{Layers}) = Layers([Layer()])

Base.:(+)(xs::Layers, ys::Layers) = Layers(vcat(xs.layers, ys.layers))

Base.:(*)(xs::Layers, y::Layer) = Layers([x*y for x in xs.layers])
Base.:(*)(xs::Layers, ys::Layers) = Layers([x*y for x in xs.layers for y in ys.layers])

function plot(args...; kw...)
    layers = one(Layers)
    for arg in args
        if arg isa Layer || arg isa Layers
            layers *= arg
        elseif arg isa AbstractVector
            layers *= Layers([plot(x) for x in arg])
        elseif arg isa Function
            layers *= Layer(; transforms=[arg])
        elseif arg isa Bokeh.ModelType && Bokeh.issubmodeltype(arg, Bokeh.Glyph)
            layers *= Layer(; glyph=arg)
        elseif arg isa Data || Bokeh.ismodelinstance(arg, Bokeh.DataSource) || istable(arg)
            layers *= plotdata(arg)
        else
            throw(ArgumentError("invalid argument of type $(typeof(arg))"))
        end
    end
    if !isempty(kw)
        layers *= Layer(; properties=Dict{Symbol,Any}(k => maybe_parse_mapping(k, v) for (k, v) in kw))
    end
    return layers
end


### DATA

function plotdata(data::Data)
    return Layers([Layer(; data)])
end

function plotdata(data::Bokeh.ModelInstance)
    if Bokeh.ismodelinstance(data, Bokeh.DataSource)
        columns = Dict{String,DataInfo}()
        if Bokeh.ismodelinstance(data, Bokeh.ColumnDataSource)
            for (name, col) in data.data
                if !isempty(col) && all(x->isa(x,Factor), col)
                    columns[name] = DataInfo(datatype=FACTOR_DATA, factors=sort(unique(col)))
                end
            end
        end
        return plotdata(Data(data, columns))
    else
        throw(ArgumentError("expecting a DataSource, got a $(Bokeh.modeltype(data).name)"))
    end
end

function plotdata(data)
    if istable(data)
        return plotdata(Bokeh.ColumnDataSource(; data))
    else
        throw(ArgumentError("expecting a table, got a $(typeof(data))"))
    end
end


### MAPPING

_recpairfirst(x) = x isa Pair ? _recpairfirst(x.first) : x
_flatpair(x) = x isa Pair ? (_flatpair(x.first)..., _flatpair(x.second)...) : (x,)

function maybe_parse_mapping(k, v)
    v0 = _recpairfirst(v)
    if v0 isa AbstractString && startswith(v0, "@")
        return parse_mapping(k, v)
    elseif v0 isa Tuple{Vararg{AbstractString}} && any(x->startswith(x, "@"), v0)
        return parse_mapping(k, v)
    elseif v0 isa Field
        return parse_mapping(k, v)
    else
        return v
    end
end

function parse_mapping(k, v)
    name = String(k)
    # infer the type from the name
    if occursin("color", name)
        type = COLOR_MAP
    elseif occursin("marker", name)
        type = MARKER_MAP
    elseif occursin("hatch_pattern", name)
        type = HATCH_PATTERN_MAP
    else
        type = DATA_MAP
    end
    # flatten any pairs into a tuple
    vs = _flatpair(v)
    # first entry is the field to map
    if vs[1] isa Field
        field = vs[1]
    elseif vs[1] isa AbstractString
        @assert startswith(vs[1], "@")
        field = Field(vs[1][2:end])
    else
        error("not implemented")
    end
    # remaining entries define optional info
    transforms = Bokeh.ModelInstance[]
    datainfo = DataInfo()
    label = nothing
    for x in vs[2:end]
        if Bokeh.ismodelinstance(x, Bokeh.Transform)
            push!(transforms, x)
        elseif x isa DataInfo
            datainfo = x
        elseif x isa AbstractString || Bokeh.ismodelinstance(x, Bokeh.BaseText)
            label = x
        else
            error("invalid mapping argument of type $(typeof(x))")
        end
    end
    return Mapping(; name, type, field, transforms, datainfo, label)
end


### DRAW

function _get_data(data, transforms, source_cache)
    if isempty(transforms)
        return data
    else
        error("not implemented")
    end
end

function _get_property(v; data)
    if v isa Mapping
        # TODO
        return v.field.name
    else
        return v
    end
end

function draw(layers::Layers)
    fig = Bokeh.figure()
    source_cache = Dict{Any,Bokeh.ModelInstance}()
    for layer in layers.layers
        # get the data
        orig_data = layer.data
        transforms = layer.transforms
        orig_data === nothing && error("no data")
        data = _get_data(orig_data, transforms, source_cache)
        source = data.source
        # get the glyph
        glyph = layer.glyph
        glyph === nothing && error("no glyph")
        # get the properties
        kw = Dict{Symbol,Any}(k => _get_property(v; data) for (k, v) in layer.properties)
        # plot the glyph
        Bokeh.plot!(fig, glyph; source, kw...)
    end
    return fig
end

function draw(; kw...)
    return layer -> draw(layer; kw...)
end

function Base.display(d::Bokeh.BokehDisplay, layers::Layers)
    # theme = get(Bokeh.setting(:theme).attrs, :Algebrokeh, nothing)
    # theme = theme === nothing ? Theme() : Theme(; theme...)
    return display(d, draw(layers))
end



# ### LAYER

# Base.@kwdef struct Layer
#     data::Any=nothing
#     transforms::Vector{Any}=[]
#     glyph::Union{Bokeh.ModelType,Nothing}=nothing
#     properties::Dict{Symbol,Any} = Dict{Symbol,Any}()
# end

# function plot(args...; kw...)
#     layers = Layers([Layer()])
#     for arg in args
#         if istable(arg) || Bokeh.ismodelinstance(arg, Bokeh.DataSource)
#             layers *= Layer(; data=arg)
#         elseif isa(arg, Bokeh.ModelType) && Bokeh.issubmodeltype(arg, Bokeh.Glyph)
#             layers *= Layer(; glyph=arg)
#         elseif isa(arg, Function)
#             layers *= Layer(; transforms=[arg])
#         elseif isa(arg, AnyLayers)
#             layers *= arg
#         elseif isa(arg, AbstractVector)
#             layers *= Layers(Layer[layer for x in arg for layer in eachlayer(plot(x))])
#         else
#             error("unexpected argument of type $(typeof(arg))")
#         end
#     end
#     if !isempty(kw)
#         layers *= Layer(; properties=kw)
#     end
#     return layers
# end

# function Base.:(*)(x::Layer, y::Layer)
#     data = y.data !== nothing ? y.data : x.data
#     glyph = y.glyph !== nothing ? y.glyph : x.glyph
#     properties = merge(x.properties, y.properties)
#     transforms = vcat(x.transforms, y.transforms)
#     return Layer(; data, glyph, properties, transforms)
# end


# ### LAYERS

# struct Layers
#     layers::Vector{Layer}
# end

# Base.convert(::Type{Layers}, x::Layers) = x
# Base.convert(::Type{Layers}, x::Layer) = Layers(Layer[x])

# Base.promote_rule(::Type{Layer}, ::Type{Layers}) = Layers

# eachlayer(x::Layer) = [x]
# eachlayer(x::Layers) = x.layers

# const AnyLayers = Union{Layer,Layers}

# function Base.:(*)(x::AnyLayers, y::AnyLayers)
#     return Layers(Layer[x*y for x in eachlayer(x) for y in eachlayer(y)])
# end

# function Base.:(+)(x::AnyLayers, y::AnyLayers)
#     return Layers(vcat(eachlayer(x), eachlayer(y)))
# end


# ### THEME

# Base.@kwdef mutable struct Theme
#     continuous_palette::Any = "Viridis"
#     categorical_palette::Any = "Dark2"
#     markers::Any = ["circle", "square", "triangle"]
#     legend_location::String = "right"
#     figure_opts::NamedTuple = NamedTuple()
# end

# function get_palette(p, n=nothing)
#     if p isa AbstractString
#         # get a named palette
#         pname = p
#         p = get(Bokeh.PALETTES, pname, nothing)
#         if p === nothing
#             # get a named palette group
#             pg = get(Bokeh.PALETTE_GROUPS, pname, nothing)
#             pg === nothing && error("no such palette: $(pname)")
#             isempty(pg) && error("empty palette group: $(pname)")
#             p = get_palette(pg, n)
#         end
#     elseif p isa AbstractDict
#         # palette group (length -> palette)
#         pg = p
#         isempty(pg) && error("empty palette group")
#         pbest = nothing
#         nbest = nothing
#         for (ncur, pcur) in pg
#             if nbest === nothing || ((n === nothing || nbest < n) ? (ncur > nbest) : (n ≤ ncur < nbest))
#                 nbest = ncur
#                 pbest = pcur
#             end
#         end
#         @assert nbest !== nothing
#         @assert pbest !== nothing
#         p = get_palette(pbest, n)
#     else
#         p = collect(String, p)
#     end
#     if n !== nothing
#         if length(p) < n
#             p = repeat(p, cld(n, length(p)))
#         end
#         if length(p) > n
#             p = p[1:n]
#         end
#         @assert length(p) == n
#     end
#     return p
# end

# function get_markers(m, n=nothing)
#     m = collect(String, m)
#     if n !== nothing
#         if length(m) < n
#             m = repeat(m, cld(n, length(m)))
#         end
#         if length(m) > n
#             m = m[1:n]
#         end
#         @assert length(m) == n
#     end
#     return m
# end


# ### MAPPING

# Base.@kwdef mutable struct Mapping
#     literal::Any = nothing
#     field::Union{String,Nothing} = nothing
#     transforms::Vector{Bokeh.ModelInstance} = Bokeh.ModelInstance[]
#     label::Any = nothing
#     datatype::Union{String,Nothing} = nothing  # number or factor
#     factors::Union{Vector{Any},Nothing} = nothing
#     mappingname::Union{String,Nothing} = nothing
#     mappingtype::Union{String,Nothing} = nothing
# end

# function is_mapping(x)
#     @nospecialize
#     if x isa AbstractString
#         return startswith(x, '@')
#     elseif x isa Pair
#         return is_mapping(x.first)
#     else
#         return x isa Bokeh.Value || x isa Bokeh.Field || x isa Bokeh.Expr
#     end
# end

# function mapfield(x::AbstractString)
#     x = convert(String, x)
#     startswith(x, '@') || error("mapping field should start with '@', got $(repr(x))")
#     return x[2:end]
# end

# function mapping!(m::Mapping, x)
#     @nospecialize
#     if x isa AbstractString
#         m.field = mapfield(x)
#         return
#     elseif x isa Bokeh.Value || x isa Bokeh.Field || x isa Bokeh.Expr
#         m.literal = x
#         return
#     elseif x isa Pair
#         if x.second isa Pair
#             mapping!(m, (x.first => x.second.first) => x.second.second)
#             return
#         elseif x.second isa AbstractString || Bokeh.ismodelinstance(x.second, Bokeh.BaseText)
#             mapping!(m, x.first)
#             m.label = x.second
#             return
#         elseif Bokeh.ismodelinstance(x.second, Bokeh.Transform)
#             mapping!(m, x.first)
#             push!(m.transforms, x.second)
#             return
#         end
#     end
#     error("invalid mapping: $(repr(x))")
# end

# function Mapping(k::Symbol, v::Any; source, theme)
#     @nospecialize
#     # get mapping type
#     mappingname = String(k)
#     if occursin("color", mappingname)
#         mappingtype = "color"
#     elseif occursin("marker", mappingname)
#         mappingtype = "marker"
#     else
#         mappingtype = "data"
#     end
#     m = Mapping(; mappingname, mappingtype)
#     # populate from the value
#     mapping!(m, v)
#     # fill in any blanks
#     if m.literal === nothing
#         m.field === nothing && error("mapping field not specified")
#         if m.datatype === nothing
#             source === nothing && error("source not specified")
#             column = source.data[m.field]
#             if all(x === missing || x isa AbstractString || x isa Tuple for x in column)
#                 m.datatype = "factor"
#             else
#                 m.datatype = "number"
#             end
#         end
#         if m.datatype == "factor" && m.factors === nothing
#             source === nothing && error("source not specified")
#             column = source.data[m.field]
#             m.factors = sort(unique(column))
#         end
#         transforms = copy(m.transforms)
#         if m.mappingtype == "data"
#             # nothing to do
#         elseif m.mappingtype == "color"
#             if m.datatype == "number"
#                 push!(transforms, Bokeh.LinearColorMapper(; palette=get_palette(theme.continuous_palette)))
#             elseif m.datatype == "factor"
#                 push!(transforms, Bokeh.CategoricalColorMapper(; palette=get_palette(theme.categorical_palette, length(m.factors)), factors=m.factors))
#             else
#                 error("unknown datatype: $(m.datatype)")
#             end
#         elseif m.mappingtype == "marker"
#             if m.datatype == "factor"
#                 push!(transforms, Bokeh.CategoricalMarkerMapper(; markers=get_markers(theme.markers, length(m.factors)), factors=m.factors))
#             else
#                 error("marker mapping must be categorical")
#             end
#         else
#             error("unknown mappingtype: $(m.datatype)")
#         end
#         if isempty(transforms)
#             m.literal = Bokeh.Field(m.field)
#         elseif length(transforms) == 1
#             m.literal = Bokeh.transform(m.field, transforms[1])
#         else
#             # TODO: https://stackoverflow.com/questions/48772907/layering-or-nesting-multiple-bokeh-transforms
#             error("not implemented: multiple transforms (on mapping $mappingname)")
#         end
#     end
#     return m
# end


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

# function linesby(cols; kw...)
#     cols = cols isa AbstractString ? String[cols] : collect(String, cols)
#     tr = let cols = cols
#         df -> combine(groupby(df, cols), [c => (c in cols ? (x -> [first(x)]) : (x -> [x])) => c for c in names(df)])
#     end
#     return plot(tr, Bokeh.MultiLine; kw...)
# end
# export linesby

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

# const MAPPING_ALIASES = Dict(
#     :x => [:xs, :right],
#     :y => [:ys, :top],
# )

# function draw!(plot::Bokeh.ModelInstance, layers::AnyLayers; theme::Theme=Theme(), legend_location="right")
#     Bokeh.ismodelinstance(plot, Bokeh.Plot) || error("expecting a Plot")
#     source_cache = IdDict{Any,Bokeh.ModelInstance}()
#     legend_cache = Dict{Any,Bokeh.ModelInstance}()
#     for layer in eachlayer(layers)
#         # get the data source
#         data = layer.data
#         transforms = layer.transforms
#         data === nothing && error("no data")
#         source = _get_source(data, transforms, source_cache)
#         # get the glyph
#         glyph = layer.glyph
#         glyph === nothing && error("no glyph")
#         # find the properties to actually use, i.e. those which are valid for the glyph
#         # and resolving aliases
#         kw = Dict{Symbol,Any}()
#         for (k, v) in layer.properties
#             if k in (:color, :alpha) || haskey(glyph.propdescs, k)
#                 kw[k] = v
#             elseif haskey(MAPPING_ALIASES, k)
#                 for k2 in MAPPING_ALIASES[k]
#                     if haskey(glyph.propdescs, k2)
#                         if !haskey(layer.properties, k2)
#                             kw[k2] = v
#                         end
#                         break
#                     end
#                 end
#             end
#         end
#         # process properties which are mappings
#         mappings = Dict{Symbol,Mapping}()
#         for (k, v) in collect(kw)
#             if is_mapping(v)
#                 m = Mapping(k, v; source, theme)
#                 mappings[k] = m
#                 kw[k] = m.literal
#             end
#         end
#         # render the glyph
#         renderer = Bokeh.plot!(plot, glyph; source, kw...)
#         # set the axis labels
#         for (k, v) in mappings
#             if k in (:x, :left, :right) && v.label !== nothing && plot.x_axis.axis_label === nothing
#                 plot.x_axis.axis_label = v.label
#             end
#             if k in (:y, :top, :bottom) && v.label !== nothing && plot.y_axis.axis_label === nothing
#                 plot.y_axis.axis_label = v.label
#             end
#         end
#         for (k, v) in mappings
#             if k in (:x, :left, :right) && v.field !== nothing && plot.x_axis.axis_label === nothing
#                 plot.x_axis.axis_label = v.field
#             end
#             if k in (:y, :top, :bottom) && v.field !== nothing && plot.y_axis.axis_label === nothing
#                 plot.y_axis.axis_label = v.field
#             end
#         end
#         # add legend items
#         for (k, v) in mappings
#             if v.datatype == "factor" && v.mappingtype in ("color", "marker")
#                 title = something(v.label, v.field, string(k))
#                 legend = get!(legend_cache, v.field) do
#                     return Bokeh.plot!(plot, Bokeh.Legend; title, location=theme.legend_location, orientation=(theme.legend_location in ("above", "below") ? "horizontal" : "vertical"))
#                 end
#                 if !any(item.label==v.literal && item.renderers==[renderer] for item in legend.items)
#                     if isempty(legend.items)
#                         push!(legend.items, Bokeh.LegendItem(label=v.literal, renderers=[renderer]))
#                     else
#                         for item in legend.items
#                             push!(item.renderers, renderer)
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return plot
# end

end # module
