using LinkedLists

include("TemporalGraph.jl")

# Cost structure
# function gamma(e)
#     return e[3]
# end
function gamma(e)
    return 1
end


# function oplus(c1, c2)
#     return c1
# end
function oplus(c1, c2)
    return c1 + c2
end

# function lessthan(c1, c2)
#     return c1 > c2
# end

# function lessthanequal(c1, c2)
#     return c1 >= c2
# end
function lessthanequal(c1, c2)
    return c1 <= c2
end

function equal(c1, c2)
    return c1 == c2
end

function lessthan(c1, c2)
    return lessthanequal(c1, c2) && !equal(c1, c2)
end

# Auxiliary functions
function distribute_edep(n, earr, edep)
    edepv = [[] for i = 1:n]
    edepv_index = zeros(Int64, length(earr))
    for i in 1:lastindex(edep)
        push!(edepv[earr[edep[i]][1]], edep[i])
        edepv_index[edep[i]] = length(edepv[earr[edep[i]][1]])
    end
    return edepv, edepv_index
end

function find_l(earr, edepv, lv, a, alphav)
    while (lv <= length(edepv))
        if (earr[edepv[lv]][3] < a + alphav)
            lv = lv + 1
        else
            break
        end
    end
    return lv
end

function find_r(earr, edepv, rv, a, betav)
    start = rv
    while (rv <= length(edepv))
        if (earr[edepv[rv]][3] > a + betav)
            break
        else
            rv = rv + 1
        end
    end
    if (rv == start)
        return rv
    else
        return rv - 1
    end
end

function print_state(e, i, a, l, r, lc, Iv, lvrv, best, sigma, cost, parent, verbose)
    if (verbose)
        logging("e: " * string(e), true, false)
        logging("i: " * string(i), true, false)
        logging("a: " * string(a), true, false)
        logging("l: " * string(l), true, false)
        logging("r: " * string(r), true, false)
        logging("l_c: " * string(lc), true, false)
        logging("I_v: " * string(Iv), true, false)
        logging("(l_v,r_v): " * string(lvrv), true, false)
        logging("B: " * string(best), true, false)
        logging("sigma: " * string(sigma), true, false)
        logging("cost: " * string(cost), true, false)
        logging("parent: " * string(parent), true, false)
        logging("====================================================", true, false)
    end
end

# Main algorithm
function finalize_cost(v, j, Iv, Pv, edepv, best, sigma, parent, lvrv)
    # println("Length of Iv: ", length(Iv))
    while (length(Iv[v]) > 0 && first(Iv[v])[1] <= j)
        I = first(Iv[v]) # I = (l,r,c,k)
        L = first(Pv[v])
        l_prime = min(I[2], j)
        for i in I[1]:l_prime
            f = edepv[v][i]
            best[f] = I[3]
            sigma[f] = I[4]
            parent[f] = L
        end
        if (l_prime == I[2])
            popfirst!(Iv[v])
            popfirst!(Pv[v])
        else
            I[1] = j + 1
        end
    end
    lvrv[v][1] = j + 1
end

function optimal_edge_walks_counter(n, alpha, beta, earr, edep, s, verbose)
    edepv, edepv_index = distribute_edep(n, earr, edep)
    logging("====================================================", verbose, false)
    logging("edepv: " * string(edepv), verbose, false)
    logging("edepv_index: " * string(edepv_index), verbose, false)
    logging("====================================================", verbose, false)
    Iv = [LinkedList{Vector{Int64}}() for i = 1:n]
    Pv = [LinkedList{Vector{Int64}}() for i = 1:n]
    lvrv = [[1, 0] for i = 1:n]
    best::Array{Union{Int64,Missing}} = [missing for e = 1:length(earr)]
    sigma = [0 for e = 1:length(earr)]
    parent = [[] for e = 1:length(earr)]
    cost::Array{Union{Int64,Missing}} = [missing for e = 1:length(earr)]
    print_state([], 0, 0, 0, 0, 0, Iv, lvrv, best, sigma, cost, parent, verbose)
    c::Int64 = 0
    for ei in 1:lastindex(earr)
        e = earr[ei]
        u = e[1]
        v = e[2]
        tau = e[3]
        lambda = e[4]
        i = edepv_index[ei]
        if (i >= lvrv[u][1])
            finalize_cost(u, i, Iv, Pv, edepv, best, sigma, parent, lvrv)
        end
        if (u == s || !ismissing(best[ei]))
            if (u == s)
                if (ismissing(best[ei]) || lessthan(gamma(e), oplus(best[ei], gamma(e))))
                    c = gamma(e)
                    sigma[ei] = 1
                end
                if (!ismissing(best[ei]) && equal(gamma(e), oplus(best[ei], gamma(e))))
                    c = gamma(e)
                    sigma[ei] = sigma[ei] + 1
                end
            else
                c = oplus(best[ei], gamma(e))
            end
            cost[ei] = c
            a::Int64 = tau + lambda
            l::Int64 = find_l(earr, edepv[v], lvrv[v][1], a, alpha[v])
            r::Int64 = find_r(earr, edepv[v], max(1, lvrv[v][2]), a, beta[v])
            finalize_cost(v, l - 1, Iv, Pv, edepv, best, sigma, parent, lvrv)
            lc = max(l, lvrv[v][2] + 1)
            while (length(Iv[v]) > 0 && lessthan(c, last(Iv[v])[3]))
                lc = last(Iv[v])[1]
                pop!(Iv[v])
            end
            if (length(Iv[v]) > 0 && equal(c, last(Iv[v])[3]))
                last(Iv[v])[4] = last(Iv[v])[4] + sigma[ei]
                push!(last(Pv[v]), ei)
            end
            if (lc <= r)
                push!(Iv[v], [lc, r, c, sigma[ei]])
                push!(Pv[v], [ei])
            end
            lvrv[v][2] = r
            if (verbose)
                print_state(e, i, a, l, r, lc, Iv, lvrv, best, sigma, cost, parent, verbose)
            end
        end
    end
    logging("====================================================", true, false)
    logging("sigma: " * string(sigma), true, false)
    logging("cost: " * string(cost), true, false)
    logging("parent: " * string(parent), true, false)
    logging("====================================================", true, false)
    return sigma, cost
end

function optimal_edge_walks_counter(fn::String, s::Int64, verbose::Bool)
    n, alpha, beta, earr, edep = read_patg(fn, ",", verbose)
    optimal_edge_walks_counter(n, alpha, beta, earr, edep, s, false)
end

# Counting number of optimal walks
# function target_cost(e, c)
#     return e[3]+e[4] - c 
# end

function target_cost(e, c)
    return c
end

function optimal_walks_counter(fn::String, s::Int64, verbose::Bool)
    n, _, _, earr, _ = read_patg(fn, ",", true)
    sigma, cost = optimal_edge_walks_counter(fn, s, false)
    optimal_value = fill(typemax(Int64), n)
    for ei in 1:lastindex(earr)
        if (!ismissing(cost[ei]))
            c = target_cost(earr[ei], cost[ei])
            if (optimal_value[earr[ei][2]] > c)
                optimal_value[earr[ei][2]] = c
            end
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging(string(optimal_value), true, false)
        logging("====================================================", true, false)
    end
    n_optimal_walks = fill(0, length(earr))
    for ei in 1:lastindex(earr)
        if (!ismissing(cost[ei]))
            c = target_cost(earr[ei], cost[ei])
            if (c == optimal_value[earr[ei][2]])
                n_optimal_walks[ei] = sigma[ei]
            end
        end
    end
    logging("====================================================", true, false)
    logging(string(n_optimal_walks), true, false)
    logging("====================================================", true, false)
    return n_optimal_walks
end
