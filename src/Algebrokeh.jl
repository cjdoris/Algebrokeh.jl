module Algebrokeh

using Base: @kwdef
using Bokeh: Bokeh
using DataFrames: DataFrame, groupby, combine, nrow
using Tables: Tables
using DataAPI: DataAPI

export plot


### TYPES

@enum DataType NUMBER_DATA FACTOR_DATA

@enum MappingType COLOR_MAP MARKER_MAP HATCH_PATTERN_MAP DATA_MAP

"""
Holds information about a collection of data, such as a column in a DataSource.
"""
@kwdef struct DataInfo
    datatype::DataType = NUMBER_DATA
    factors::Vector = []
    label::Any = nothing
end

"""
A DataSource, plus information about the data in its columns.
"""
@kwdef struct Data
    source::Union{Bokeh.ModelInstance}
    table::Any
    columns::Dict{String,DataInfo}
end

"""
A specifier for a field. Can also be a collection of 2 or 3 fields, to treat them as a
hierarchical factor.
"""
struct Field
    names::Union{String,Tuple{String,String},Tuple{String,String,String}}
end

"""
A mapping value, which defines some connection between a field of data and a visual
property, such as color, marker or location.
"""
@kwdef struct Mapping
    name::String
    type::MappingType
    field::Field
    label::Any = nothing
    transforms::Vector{Bokeh.ModelInstance} = Bokeh.ModelInstance[]
    datainfo::Union{DataInfo,Nothing} = nothing
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

"""
    plot(args...; kw...)

Create an Algebrokeh plot.
"""
function plot(args...; kw...)
    layers = one(Layers)
    for arg in args
        if arg isa Layer || arg isa Layers
            layers *= arg
        elseif arg isa AbstractVector
            layers *= Layers([layer for x in arg for layer in plot(x).layers])
        elseif arg isa Function
            layers *= Layer(; transforms=[arg])
        elseif arg isa Bokeh.ModelType && Bokeh.issubmodeltype(arg, Bokeh.Glyph)
            layers *= Layer(; glyph=arg)
        elseif arg isa Data || Bokeh.ismodelinstance(arg, Bokeh.DataSource) || Tables.istable(arg)
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

function plotdata(data)
    if !isa(data, Data)
        if Bokeh.ismodelinstance(data, Bokeh.ColumnDataSource)
            table = DataFrame(data.data)
            source = data
        elseif Bokeh.ismodelinstance(data, Bokeh.DataSource)
            table = nothing
            source = data
        elseif Bokeh.ismodelinstance(data)
            throw(ArgumentError("expecting a DataSource, got a $(Bokeh.modeltype(data).name)"))
        elseif Tables.istable(data)
            table = data
            source = Bokeh.ColumnDataSource(; data)
        else
            throw(ArgumentError("expecting a table, got a $(typeof(data))"))
        end
        columns = Dict{String,DataInfo}()
        if table !== nothing
            tablecols = Tables.columns(table)
            for name in Tables.columnnames(tablecols)
                col = Tables.getcolumn(tablecols, name)
                if DataAPI.colmetadatasupport(typeof(table)).read
                    label = DataAPI.colmetadata(table, name, "label", nothing)
                else
                    label = nothing
                end
                if !isempty(col) && all(x->isa(x, Factor), col)
                    datatype = FACTOR_DATA
                    factors = sort(unique(col))
                else
                    datatype = NUMBER_DATA
                    factors = []
                end
                columns[String(name)] = DataInfo(; datatype, factors, label)
            end
        end
        data = Data(; source, table, columns)
    end
    return Layer(; data)
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
    datainfo = nothing
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

function _get_palette(p, n=nothing)
    if p isa AbstractString
        # get a named palette
        pname = p
        p = get(Bokeh.PALETTES, pname, nothing)
        if p === nothing
            # get a named palette group
            pg = get(Bokeh.PALETTE_GROUPS, pname, nothing)
            pg === nothing && error("no such palette: $(pname)")
            isempty(pg) && error("empty palette group: $(pname)")
            p = _get_palette(pg, n)
        end
    elseif p isa AbstractDict
        # palette group (length -> palette)
        pg = p
        isempty(pg) && error("empty palette group")
        pbest = nothing
        nbest = nothing
        for (ncur, pcur) in pg
            if nbest === nothing || ((n === nothing || nbest < n) ? (ncur > nbest) : (n ≤ ncur < nbest))
                nbest = ncur
                pbest = pcur
            end
        end
        @assert nbest !== nothing
        @assert pbest !== nothing
        p = _get_palette(pbest, n)
    else
        p = collect(String, p)
    end
    if n !== nothing
        if length(p) < n
            p = repeat(p, cld(n, length(p)))
        end
        if length(p) > n
            p = p[1:n]
        end
        @assert length(p) == n
    end
    return p
end

function _get_markers(m, n=nothing)
    m = collect(String, m)
    if n !== nothing
        if length(m) < n
            m = repeat(m, cld(n, length(m)))
        end
        if length(m) > n
            m = m[1:n]
        end
        @assert length(m) == n
    end
    return m
end

function _get_hatch_patterns(p, n=nothing)
    return _get_markers(p, n)
end

function _get_data(data, transforms, source_cache)
    if isempty(transforms)
        return data
    else
        error("not implemented")
    end
end

@kwdef struct ResolvedProperty
    orig::Any
    value::Any
    label::Any = nothing
    field::Union{Field,Nothing} = nothing
    fieldname::Union{String,Nothing} = nothing
    datainfo::Union{DataInfo,Nothing} = nothing
end

function _get_property(v; data::Data)
    if v isa Mapping
        # field
        field = v.field
        fieldnames = field.names
        if fieldnames isa String
            fieldname = fieldnames
        else
            # TODO: add a new combined column
            error("not implemented: hierarchical fields")
        end
        # datainfo
        datainfo = v.datainfo
        if datainfo === nothing
            datainfo = get(data.columns, fieldname, nothing)
            if datainfo === nothing
                error("no DataInfo for column $(repr(field))")
            end
        end
        datatype = datainfo.datatype
        factors = datainfo.factors
        nfactors = length(factors)
        label = v.label
        if label === nothing
            label = datainfo.label
        end
        # transforms
        transforms = copy(v.transforms)
        if v.type == DATA_MAP
            # no extra transforms
        elseif v.type == COLOR_MAP
            if datatype == NUMBER_DATA
                palette = _get_palette("Viridis256")
                push!(transforms, Bokeh.LinearColorMapper(; palette))
            else @assert datatype == FACTOR_DATA
                palette = _get_palette("Dark2", nfactors)
                push!(transforms, Bokeh.CategoricalColorMapper(; palette, factors))
            end
        elseif v.type == MARKER_MAP
            if datatype == FACTOR_DATA
                markers = _get_markers(["circle", "square", "triangle"], nfactors)
                push!(transforms, Bokeh.CategoricalMarkerMapper(; markers, factors))
            else
                error("$(v.name) is a marker mapping but $fields is not categorical")
            end
        else @assert v.type == HATCH_PATTERN_MAP
            if datatype == FACTOR_DATA
                patterns = _get_hatch_patterns(["/", "\\", "+", ".", "o"], nfactors)
                push!(transforms, Bokeh.CategoricalPatternMapper(; patterns, factors))
            else
                error("$(v.name) is a hatch-pattern mapping but $fields is not categorical")
            end
        end
        # done
        if length(transforms) == 0
            value = Bokeh.Field(fieldname)
        elseif length(transforms) == 1
            value = Bokeh.transform(fieldname, transforms[1])
        else
            # TODO: https://stackoverflow.com/questions/48772907/layering-or-nesting-multiple-bokeh-transforms
            error("not implemented: multiple transforms (on mapping $(v.name))")
        end
        return ResolvedProperty(; orig=v, value, label, field, fieldname, datainfo)
    else
        return ResolvedProperty(; orig=v, value=v)
    end
end

const PROPERTY_ALIASES = Dict(
    :x => [:xs, :right],
    :y => [:ys, :top],
)

function _get_label(layers, keys)
    props = Any[v for layer in layers for (k, v) in layer.props if k in keys]
    labels = Any[p.label for p in props if p.label !== nothing]
    if !isempty(labels)
        return first(labels)
    end
    labels = String[string(p.field.names) for p in props if p.field !== nothing]
    if !isempty(labels)
        return join(sort(unique(labels)), " / ")
    end
    return nothing
end

@kwdef struct ResolvedLayer
    orig::Layer
    data::Data
    props::Dict{Symbol,ResolvedProperty}
    renderer::Bokeh.ModelInstance
end

function draw(layers::Layers)
    fig = Bokeh.figure()

    # PLOT/RESOLVE EACH LAYER
    source_cache = Dict{Any,Bokeh.ModelInstance}()
    resolved = Vector{ResolvedLayer}()
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
        # get the properties actually applicable to this type, resolving aliases too
        props = Dict{Symbol,Any}()
        for (k, v) in layer.properties
            if k in (:color, :alpha) || haskey(glyph.propdescs, k)
                props[k] = v
            elseif haskey(PROPERTY_ALIASES, k)
                for k2 in PROPERTY_ALIASES[k]
                    if haskey(glyph.propdescs, k2)
                        if !haskey(layer.properties, k2)
                            props[k2] = v
                        end
                        break
                    end
                end
            end
        end
        # resolve the properties
        props = Dict{Symbol,ResolvedProperty}(k => _get_property(v; data) for (k, v) in props)
        # plot the glyph
        kw = Dict{Symbol,Any}(k => v.value for (k, v) in props)
        renderer = Bokeh.plot!(fig, glyph; source, kw...)
        # save
        push!(resolved, ResolvedLayer(; orig=layer, data, props, renderer))
    end

    # AXIS LABELS
    fig.x_axis.axis_label = _get_label(resolved, [:x, :xs, :right, :left])
    fig.y_axis.axis_label = _get_label(resolved, [:y, :ys, :top, :bottom])

    # LEGENDS
    # Find all layers+props for categorical mappings.
    legendinfos = Dict()
    for layer in resolved
        for (k, v) in layer.props
            v.datainfo !== nothing || continue
            v.datainfo.datatype == FACTOR_DATA || continue
            m = v.orig
            m isa Mapping || continue
            m.type in (COLOR_MAP, MARKER_MAP, HATCH_PATTERN_MAP) || continue
            push!(get!(legendinfos, (layer.data.source, v.fieldname), []), (layer, v))
        end
    end
    # Generate a legend for each unique source+field combination.
    legends = []
    for ((source, fieldname), layerprops) in legendinfos
        item = Bokeh.LegendItem(; label=Bokeh.Field(fieldname), renderers=[layer.renderer for (layer, _) in layerprops])
        items = [item]
        titles = [p.label for (_, p) in layerprops if p.label !== nothing]
        if !isempty(titles)
            title = titles[1]
        else
            title = fieldname
        end
        legend = Bokeh.Legend(; items, title)
        push!(legends, (fieldname, legend))
    end
    # Plot the legends, sorted by title.
    sort!(legends, by=x->x[1])
    for (_, legend) in legends
        Bokeh.plot!(fig, legend, location="right")
    end

    # TODO: COLOR BARS

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
