function linesby(cols...; kw...)
    cols = collect(String, cols)
    function tr(data::Data)
        table0 = data.table
        table0 === nothing && error("linesby() requires data to be specified explicitly")
        df = DataFrame(table0)
        table = combine(groupby(df, cols), [c => (c in cols ? (x -> [first(x)]) : (x -> [x])) => c for c in names(df)])
        source = Bokeh.ColumnDataSource(; data=table)
        columns = copy(data.columns)
        return Data(; table, source, columns)
    end
    return plot(tr, Bokeh.MultiLine; kw...)
end
