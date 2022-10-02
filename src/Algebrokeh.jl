module Algebrokeh

import Bokeh

export data, visual, mapping, vstack, hstack, draw!, draw


### LAYER

Base.@kwdef struct Layer
    data::Any=nothing
    glyph::Union{Bokeh.ModelType,Nothing}=nothing
    mappings::Dict{Symbol,Any} = Dict{Symbol,Any}()
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

data(data) = Layer(; data)
visual(glyph::Bokeh.ModelType; kw...) = Layer(; glyph, properties=Dict{Symbol,Any}(kw))
visual(; kw...) = Layer(; properties=Dict{Symbol,Any}(kw))
mapping(; kw...) = Layer(; mappings=Dict{Symbol,Any}(kw))

function Base.:(*)(x::Layer, y::Layer)
    data = y.data !== nothing ? y.data : x.data
    glyph = y.glyph !== nothing ? y.glyph : x.glyph
    mappings = merge(x.mappings, y.mappings)
    properties = merge(x.properties, y.properties)
    return Layer(; data, glyph, mappings, properties)
end


### LAYERS

struct Layers
    layers::Vector{Layer}
end

Base.convert(::Type{Layers}, x::Layers) = x
Base.convert(::Type{Layers}, x::Layer) = Layers(Layer[x])

eachlayer(x::Layer) = [x]
eachlayer(x::Layers) = x.layers

const AnyLayers = Union{Layer,Layers}

function Base.:(*)(x::AnyLayers, y::AnyLayers)
    return Layers(Layer[x*y for x in eachlayer(x) for y in eachlayer(y)])
end

function Base.:(+)(x::AnyLayers, y::AnyLayers)
    return Layers(vcat(eachlayer(x), eachlayer(y)))
end


### THEME

Base.@kwdef mutable struct Theme
    continuous_palette::Any = "Viridis"
    categorical_palette::Any = "Dark2"
    markers::Any = ["circle", "square", "triangle"]
end

function get_palette(p, n=nothing)
    if p isa AbstractString
        # get a named palette
        pname = p
        p = get(Bokeh.PALETTES, pname, nothing)
        if p === nothing
            # get a named palette group
            pg = get(Bokeh.PALETTE_GROUPS, pname, nothing)
            pg === nothing && error("no such palette: $(pname)")
            isempty(pg) && error("empty palette group: $(pname)")
            p = get_palette(pg, n)
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
        p = get_palette(pbest, n)
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

function get_markers(m, n=nothing)
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


### MAPPING

Base.@kwdef mutable struct Mapping
    literal::Any = nothing
    field::Union{String,Nothing} = nothing
    label::Any = nothing
    datatype::Union{String,Nothing} = nothing  # number or factor
    factors::Union{Vector{Any},Nothing} = nothing
    mappingname::Union{String,Nothing} = nothing
    mappingtype::Union{String,Nothing} = nothing
end

function mapping!(m::Mapping, x)
    @nospecialize
    if x isa AbstractString
        m.field = x
        return
    elseif x isa Bokeh.Value || x isa Bokeh.Field || x isa Bokeh.Expr
        m.literal = x
        return
    elseif x isa Pair
        if x.second isa Pair
            mapping!(m, (x.first => x.second.first) => x.second.second)
            return
        elseif x.second isa AbstractString || Bokeh.ismodelinstance(x.second, Bokeh.BaseText)
            mapping!(m, x.first)
            m.label = x.second
            return
        end
    end
    error("invalid mapping: $(repr(x))")
end

function Mapping(k::Symbol, v::Any; source, theme)
    @nospecialize
    # get mapping type
    mappingname = String(k)
    if occursin("color", mappingname)
        mappingtype = "color"
    elseif occursin("marker", mappingname)
        mappingtype = "marker"
    else
        mappingtype = "data"
    end
    m = Mapping(; mappingname, mappingtype)
    # populate from the value
    mapping!(m, v)
    # fill in any blanks
    if m.literal === nothing
        m.field === nothing && error("mapping field not specified")
        if m.datatype === nothing
            source === nothing && error("source not specified")
            column = source.data[m.field]
            if all(x === missing || x isa AbstractString || x isa Tuple for x in column)
                m.datatype = "factor"
            else
                m.datatype = "number"
            end
        end
        if m.datatype == "factor" && m.factors === nothing
            source === nothing && error("source not specified")
            column = source.data[m.field]
            m.factors = sort(unique(column))
        end
        if m.mappingtype == "data"
            m.literal = Bokeh.Field(m.field)
        elseif m.mappingtype == "color"
            if m.datatype == "number"
                m.literal = Bokeh.linear_cmap(m.field, get_palette(theme.continuous_palette))
            elseif m.datatype == "factor"
                m.literal = Bokeh.factor_cmap(m.field, get_palette(theme.categorical_palette, length(m.factors)), m.factors)
            else
                error("unknown datatype: $(m.datatype)")
            end
        elseif m.mappingtype == "marker"
            if m.datatype == "factor"
                m.literal = Bokeh.factor_mark(m.field, get_markers(theme.markers, length(m.factors)), m.factors)
            else
                error("marker mapping must be categorical")
            end
        else
            error("unknown mappingtype: $(m.datatype)")
        end
    end
    return m
end


# ### VSTACK

function stack(k1, k2, fields; kw...)
    fields = collect(String, fields)
    layers = Layer[]
    x1 = Bokeh.Expr(Bokeh.Stack(fields=[]))
    for i in 1:length(fields)
        x2 = Bokeh.Expr(Bokeh.Stack(fields=fields[1:i]))
        kwnew = Dict{Symbol,Any}()
        for (k, v) in kw
            if v isa AbstractVector
                kwnew[k] = v[i]
            else
                kwnew[k] = v
            end
        end
        layer = Layer(mappings=Dict(k1=>x1, k2=>x2), properties=kwnew)
        push!(layers, layer)
        x1 = x2
    end
    return Layers(layers)
end

function vstack(fields; kw...)
    return stack(:bottom, :top, fields; kw...)
end

function hstack(fields; kw...)
    return stack(:left, :right, fields; kw...)
end


### DRAW

function draw!(plot::Bokeh.ModelInstance, layers::AnyLayers; theme::Theme=Theme())
    Bokeh.ismodelinstance(plot, Bokeh.Plot) || error("expecting a Plot")
    source_cache = IdDict()
    legend_cache = Dict()
    for layer in eachlayer(layers)
        # get the data source
        data = layer.data
        data === nothing && error("no data")
        source = get!(source_cache, data) do
            return Bokeh.ColumnDataSource(; data)
        end
        # get the glyph
        glyph = layer.glyph
        glyph === nothing && error("no glyph")
        # get the mappings actually used
        orig_mappings = Dict(k => Mapping(k, v; source, theme) for (k, v) in layer.mappings)
        mappings = empty(orig_mappings)
        for (k, v) in orig_mappings
            if k in (:color, :alpha) || haskey(glyph.propdescs, k)
                mappings[k] = v
            elseif k == :x && !haskey(orig_mappings, :xs) && haskey(glyph.propdescs, :xs)
                mappings[:xs] = v
            elseif k == :y && !haskey(orig_mappings, :ys) && haskey(glyph.propdescs, :ys)
                mappings[:ys] = v
            end
        end
        # collect all glyph properties and render the glyph
        kw = collect(layer.properties)
        for (k, v) in mappings
            push!(kw, Pair{Symbol,Any}(k, v.literal))
        end
        renderer = Bokeh.plot!(plot, glyph; source, kw...)
        # set the axis labels
        for (k, v) in mappings
            if k in (:x, :left, :right) && v.label !== nothing && plot.x_axis.axis_label === nothing
                plot.x_axis.axis_label = v.label
            end
            if k in (:y, :top, :bottom) && v.label !== nothing && plot.y_axis.axis_label === nothing
                plot.y_axis.axis_label = v.label
            end
        end
        for (k, v) in mappings
            if k in (:x, :left, :right) && v.field !== nothing && plot.x_axis.axis_label === nothing
                plot.x_axis.axis_label = v.field
            end
            if k in (:y, :top, :bottom) && v.field !== nothing && plot.y_axis.axis_label === nothing
                plot.y_axis.axis_label = v.field
            end
        end
        # add legend items
        for (k, v) in mappings
            if v.datatype == "factor" && v.mappingtype in ("color", "marker")
                title = something(v.label, v.field, string(k))
                legend = get!(legend_cache, title) do
                    return Bokeh.plot!(plot, Bokeh.Legend; title, location="right")
                end
                if !any(item.label==v.literal && item.renderers==[renderer] for item in legend.items)
                    if isempty(legend.items)
                        push!(legend.items, Bokeh.LegendItem(label=v.literal, renderers=[renderer]))
                    else
                        for item in legend.items
                            push!(item.renderers, renderer)
                        end
                    end
                end
            end
        end
    end
    return plot
end

function draw(layer::AnyLayers; figure=NamedTuple(), theme=NamedTuple())
    figure = figure isa Bokeh.ModelInstance ? figure : Bokeh.figure(; figure...)
    theme = theme isa Theme ? theme : Theme(; theme...)
    draw!(figure, layer; theme)
    return figure
end

end # module
