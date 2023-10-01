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

function read_patg(fn::String, sep::String, verbose::Bool)
    @assert isfile(fn) "The point availability temporal graph file " * fn * " does not exist"
    f::IOStream = open(fn, "r")
    n::Int64 = parse(Int64, readline(f))
    alpha::Array{Int64} = fill(0, n)
    beta::Array{Int64} = fill(typemax(Int64) ÷ 2, n)
    for v in 1:n
        l::String = readline(f)
        sl::Array{String} = split(l, sep)
        @assert length(sl) >= 1 "Bad line format: " * l
        alpha[v] = parse(Int64, sl[1])
        if (length(sl) > 1)
            beta[v] = parse(Int64, sl[2])
        end
        @assert beta[v] >= alpha[v] "Beta value smaller than alpha value: " * v
    end
    earr::Vector{Array{Int64}} = []
    while (!eof(f))
        l::String = readline(f)
        sl::Array{String} = split(l, sep)
        @assert length(sl) == 4 "Bad line format: " * l
        u::Int64 = parse(Int64, sl[1])
        v::Int64 = parse(Int64, sl[2])
        tau::Int64 = parse(Int64, sl[3])
        lambda::Int64 = parse(Int64, sl[4])
        push!(earr, [u, v, tau, lambda])
    end
    sort!(earr, by=e -> e[3] + e[4])
    edep = sortperm(earr, by=e -> e[3])
    if (verbose)
        print_patg(n, alpha, beta, earr, edep)
    end
    return n, alpha, beta, earr, edep
end

function print_patg(n, alpha, beta, earr, edep)
    logging("====================================================", true, false)
    logging("Point availability temporal graph", true, false)
    logging("====================================================", true, false)
    logging("Number of nodes: " * string(n), true, false)
    logging("Temporal edges sorted by arrival time: " * string(earr), true, false)
    logging("Temporal edges sorted by departure time: " * string(edep), true, false)
    logging("Waiting lower bound constraints: " * string(alpha), true, false)
    logging("Waiting upper bound constraints: " * string(beta), true, false)
    logging("====================================================", true, false)

end
