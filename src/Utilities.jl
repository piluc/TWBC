# COMMON FUNCTIONS
function logging(s::String)
    println(s)
    flush(stdout)
end

function logging(s::String, verbose::Bool, time::Bool)
    if (verbose)
        if (time)
            println(now(), " ", s)
        else
            println(s)
        end
        flush(stdout)
    end
end

function read_centrality_values(file_name::String)::Array{Float64}
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = []
    value::Float64 = 0.0
    for l in eachline(f)
        value = parse(Float64, l)
        if (value < -0.1)
            logging("ERROR. There are negative values with absolute big values (" * l * ")")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, value)
    end
    close(f)
    return centrality
end

function save_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "w")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end

function compare_centrality_values(b1, b2, abs_tol, rel_tol)
    for i in 1:lastindex(b1)
        if (!isapprox(b1[i], b2[i]; atol=abs_tol, rtol=rel_tol))
            logging(string(b1[i]) * " != " * string(b2[i]) * " at " * string(i), true, false)
            return
        end
    end
    logging("Centrality values are the same", true, false)
end

# DATA STRUCTURE AND FUNCTIONS FOR BCV
struct PATG
    n::Int64
    α::Array{Int64}
    β::Array{Int64}
    earr::Vector{Array{Int64}}
    edep::Vector{Int64}
    edepv::Array{Vector{Int64}}
    edepv_index::Array{Int64}
    edepv_dep::Array{Vector{Int64}}
    function PATG(_n, _α, _β, _earr, _edep, _edepv, _edepv_index, _edepv_dep)
        return new(_n, _α, _β, _earr, _edep, _edepv, _edepv_index, _edepv_dep)
    end
end

function distribute_edep(n, earr, edep)
    edepv = [[] for i = 1:n]
    edepv_index = zeros(Int64, length(earr))
    for i in 1:lastindex(edep)
        push!(edepv[earr[edep[i]][1]], edep[i])
        edepv_index[edep[i]] = length(edepv[earr[edep[i]][1]])
    end
    return edepv, edepv_index
end

function departures(n, earr, edepv)
    edepv_dep::Array{Vector{Int64}} = [[] for i = 1:n]
    for v in 1:n
        for ei in edepv[v]
            push!(edepv_dep[v], earr[ei][3])
        end
    end
    return edepv_dep
end

function read_patg(fn::String, sep::String; α=0, β=typemax(Int64) ÷ 2)
    @assert isfile(fn) "The point availability temporal graph file " * fn * " does not exist"
    f::IOStream = open(fn, "r")
    l::String = readline(f)
    sl::Array{String} = split(l, sep)
    @assert length(sl) == 1 "Bad line format: " * l
    n::Int64 = parse(Int64, sl[1])
    α_v::Array{Int64} = fill(α, n)
    β_v::Array{Int64} = fill(β, n)
    if (α >= 0)
        @assert β >= α "Beta value smaller than alpha value: " * string(β) * "<" * string(α)
    else
        for v in 1:n
            l = readline(f)
            sl = split(l, sep)
            @assert length(sl) >= 1 && length(sl) <= 2 "Bad line format: " * l
            α_v[v] = parse(Int64, sl[1])
            if (length(sl) > 1)
                β_v[v] = parse(Int64, sl[2])
            else
                β_v[v] = typemax(Int64) ÷ 2
            end
            @assert β_v[v] >= α_v[v] "Beta value smaller than alpha value: " * string(v)
        end
    end
    earr::Vector{Array{Int64}} = []
    while (!eof(f))
        l = readline(f)
        sl = split(l, sep)
        @assert (length(sl) >= 3 && length(sl) < 5) "Bad line format: " * l * " " * string(sl)
        u::Int64 = parse(Int64, sl[1])
        v::Int64 = parse(Int64, sl[2])
        τ::Int64 = parse(Int64, sl[3])
        λ::Int64 = 1
        if (length(sl) == 4)
            λ = parse(Int64, sl[4])
        end
        push!(earr, [u, v, τ, λ])
    end
    close(f)
    sort!(earr, by=e -> e[3] + e[4])
    edep = sortperm(earr, by=e -> e[3])
    edepv, edepv_index = distribute_edep(n, earr, edep)
    edepv_dep = departures(n, earr, edepv)
    return PATG(n, α_v, β_v, earr, edep, edepv, edepv_index, edepv_dep)
end

function print_patg_stats(tg)
    duv = Dict{Array{Int64},Bool}()
    dt = Dict{Int64,Bool}()
    for e in tg.earr
        uv = [e[1], e[2]]
        duv[uv] = true
        dt[e[3]] = true
    end
    logging("====================================================", true, false)
    logging("Point availability temporal graph", true, false)
    logging("====================================================", true, false)
    logging("Number of nodes: " * string(tg.n), true, false)
    logging("Number of edges: " * string(length(duv)), true, false)
    logging("Number of temporal edges: " * string(length(tg.earr)), true, false)
    logging("Number of distinct time steps: " * string(length(dt)), true, false)
    logging("====================================================", true, false)
end

# COST STRUCTURES FOR BCV (ALGORITHM 4)
struct DH_COST
    dep::Int64
    hop::Int64
    function DH_COST(_d, _h)
        return new(_d, _h)
    end
end

@enum CS TSB = 1 TFaB = 2 TFoB = 3 TSFoB = 4 TSFaB = 5 TLaB = 6 TSLaB = 7

function TSB_CT()
    return TSB
end

function TFaB_CT()
    return TFaB
end

function TFoB_CT()
    return TFoB
end

function TSFoB_CT()
    return TSFoB
end

function TSFaB_CT()
    return TSFaB
end

function TLaB_CT()
    return TLaB
end

function TSLaB_CT()
    return TSLaB
end

function cs_of_string(c::String)
    if c == "S"
        return TSB
    end
    if c == "Fa"
        return TFaB
    end
    if c == "Fo"
        return TFoB
    end
    if c == "SFo"
        return TSFoB
    end
    if c == "SFa"
        return TSFaB
    end
    if c == "L"
        return TLaB
    end
    if c == "SL"
        return TSLaB
    end
    @assert false "Unknown cost"
end

@inline function cs_edge_cost(ct, e)
    return DH_COST(e[3], 1)
end

@inline function cs_oplus(ct, c1, c2)
    return DH_COST(c1.dep, c1.hop + c2.hop)
end

@inline function cs_lessthan(ct, c1, c2)
    if ct == TSB
        return c1.hop < c2.hop
    end
    if ct == TFaB
        return c1.dep > c2.dep
    end
    if ct == TFoB
        return false
    end
    if ct == TSFoB
        return c1.hop < c2.hop
    end
    if ct == TLaB
        return false
    end
    if ct == TSLaB
        return c1.hop < c2.hop
    end
    if ct == TSFaB
        return c1.dep > c2.dep || (c1.dep == c2.dep && c1.hop < c2.hop)
    end
    @assert false "Unknown cost"
end

@inline function cs_equal(ct, c1, c2)
    if ct == TSB
        return c1.hop == c2.hop
    end
    if ct == TFaB
        return c1.dep == c2.dep
    end
    if ct == TFoB
        return true
    end
    if ct == TSFoB
        return c1.hop == c2.hop
    end
    if ct == TLaB
        return true
    end
    if ct == TSLaB
        return c1.hop == c2.hop
    end
    if ct == TSFaB
        return c1.dep == c2.dep && c1.hop == c2.hop
    end
    @assert false "Unknown cost"
end

@inline function cs_lessthanequal(ct, c1, c2)
    if ct == TSB
        return c1.hop <= c2.hop
    end
    if ct == TFaB
        return c1.dep >= c2.dep
    end
    if ct == TFoB
        return true
    end
    if ct == TSFoB
        return c1.hop <= c2.hop
    end
    if ct == TLaB
        return true
    end
    if ct == TSLaB
        return c1.hop <= c2.hop
    end
    if ct == TSFaB
        return c1.dep > c2.dep || (c1.dep == c2.dep && c1.hop <= c2.hop)
    end
    @assert false "Unknown cost"
end

struct TARGET_COST
    dur::Int64
    arr::Int64
    hop::Int64
    function TARGET_COST(_d, _a, _h)
        return new(_d, _a, _h)
    end
end

@inline function cs_min_target_cost(ct)
    if ct == TLaB || ct == TSLaB
        TARGET_COST(0, typemax(Int64), 0)
    else
        TARGET_COST(0, typemin(Int64), 0)
    end
end

@inline function cs_max_target_cost(ct)
    if ct == TLaB || ct == TSLaB
        TARGET_COST(typemax(Int64), typemin(Int64), typemax(Int64))
    else
        TARGET_COST(typemax(Int64), typemax(Int64), typemax(Int64))
    end
end

@inline function cs_target_cost(ct, e, c)
    arr = e[3] + e[4]
    return TARGET_COST(arr - c.dep, arr, c.hop)
end

@inline function cs_target_short(ct, c)
    if ct == TSB
        return c.hop
    end
    if ct == TFaB
        return c.dur
    end
    if ct == TFoB
        return c.arr
    end
    if ct == TSFoB
        return (c.arr, c.hop)
    end
    if ct == TLaB
        return c.arr
    end
    if ct == TSLaB
        return (c.arr, c.hop)
    end
    if ct == TSFaB
        return (c.dur, c.hop)
    end
    @assert false "Unknown cost"
end

@inline function cs_target_infty(ct)
    if ct == TSB
        return typemax(Int64)
    end
    if ct == TFaB
        return typemax(Int64)
    end
    if ct == TFoB
        return typemax(Int64)
    end
    if ct == TSFoB
        return (typemax(Int64), typemax(Int64))
    end
    if ct == TLaB
        return typemin(Int64)
    end
    if ct == TSLaB
        return (typemin(Int64), typemax(Int64))
    end
    if ct == TSFaB
        return (typemax(Int64), typemax(Int64))
    end
    @assert false "Unknown cost"
end

@inline function cs_target_lessthan(ct, c1, c2)
    if ct == TSB
        return c1.hop < c2.hop
    end
    if ct == TFaB
        return c1.dur < c2.dur
    end
    if ct == TFoB
        return c1.arr < c2.arr
    end
    if ct == TSFoB
        return c1.arr < c2.arr || (c1.arr == c2.arr && c1.hop < c2.hop)
    end
    if ct == TLaB
        return c1.arr > c2.arr
    end
    if ct == TSLaB
        return c1.arr > c2.arr || (c1.arr == c2.arr && c1.hop < c2.hop)
    end
    if ct == TSFaB
        return c1.dur < c2.dur || (c1.dur == c2.dur && c1.hop < c2.hop)
    end
    @assert false "Unknown cost"
end

@inline function cs_target_equal(ct, c1, c2)
    if ct == TSB
        return c1.hop == c2.hop
    end
    if ct == TFaB
        return c1.dur == c2.dur
    end
    if ct == TFoB
        return c1.arr == c2.arr
    end
    if ct == TSFoB
        return c1.arr == c2.arr && c1.hop == c2.hop
    end
    if ct == TLaB
        return c1.arr == c2.arr
    end
    if ct == TSLaB
        return c1.arr == c2.arr && c1.hop == c2.hop
    end
    if ct == TSFaB
        return c1.dur == c2.dur && c1.hop == c2.hop
    end
    @assert false "Unknown cost"
end

