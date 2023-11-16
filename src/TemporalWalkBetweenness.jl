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
function distribute_edep(n, earr, edep, verbose)
    edepv = [[] for i = 1:n]
    edepv_index = zeros(Int64, length(earr))
    for i in 1:lastindex(edep)
        push!(edepv[earr[edep[i]][1]], edep[i])
        edepv_index[edep[i]] = length(edepv[earr[edep[i]][1]])
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("edepv: " * string(edepv), true, false)
        logging("edepv_index: " * string(edepv_index), true, false)
        logging("====================================================", true, false)
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
    cur = max(1, rv)
    one_found = false
    while (cur <= length(edepv))
        if (earr[edepv[cur]][3] - a > betav)
            break
        else
            one_found = true
            cur = cur + 1
        end
    end
    if (!one_found)
        return rv
    else
        return cur - 1
    end
end

function print_state(e, i, a, l, r, lc, edepv, edepv_index, Iv, Pv, lvrv, best, sigma, cost, parent, verbose)
    if (verbose)
        logging("====================================================", true, false)
        logging("e: " * string(e), true, false)
        logging("i: " * string(i), true, false)
        logging("a: " * string(a), true, false)
        logging("l: " * string(l), true, false)
        logging("r: " * string(r), true, false)
        logging("l_c: " * string(lc), true, false)
        logging("edepv: " * string(edepv), true, false)
        logging("edepv_index: " * string(edepv_index), true, false)
        logging("I_v: " * string(Iv), true, false)
        logging("P_v: " * string(Pv), true, false)
        logging("(l_v,r_v): " * string(lvrv), true, false)
        logging("best: " * string(best), true, false)
        logging("sigma: " * string(sigma), true, false)
        logging("cost: " * string(cost), true, false)
        logging("parent: " * string(parent), true, false)
        logging("====================================================", true, false)
    end
end

# Counting number of optimal walks ending with edge
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

function optimal_walks_counter(n, alpha, beta, earr, edepv, edepv_index, s, verbose)
    Iv = [LinkedList{Vector{Int64}}() for i = 1:n]
    Pv = [LinkedList{Vector{Int64}}() for i = 1:n]
    lvrv = [[1, 0] for i = 1:n]
    best::Array{Union{Int64,Missing}} = [missing for e = 1:length(earr)]
    sigma = [0 for e = 1:length(earr)]
    parent = [[] for e = 1:length(earr)]
    cost::Array{Union{Int64,Missing}} = [missing for e = 1:length(earr)]
    print_state([], 0, 0, 0, 0, 0, edepv, edepv_index, Iv, Pv, lvrv, best, sigma, cost, parent, verbose)
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
                    parent[ei] = []
                end
                if (!ismissing(best[ei]) && equal(gamma(e), oplus(best[ei], gamma(e))))
                    c = gamma(e)
                    sigma[ei] = sigma[ei] + 1
                    # What to do? Still a problem!
                end
            else
                c = oplus(best[ei], gamma(e))
            end
            cost[ei] = c
            a::Int64 = tau + lambda
            l::Int64 = find_l(earr, edepv[v], lvrv[v][1], a, alpha[v])
            r::Int64 = find_r(earr, edepv[v], lvrv[v][2], a, beta[v])
            finalize_cost(v, l - 1, Iv, Pv, edepv, best, sigma, parent, lvrv)
            lc = max(l, lvrv[v][2] + 1)
            while (length(Iv[v]) > 0 && lessthan(c, last(Iv[v])[3]))
                lc = last(Iv[v])[1]
                pop!(Iv[v])
                pop!(Pv[v])
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
                print_state(e, i, a, l, r, lc, edepv, edepv_index, Iv, Pv, lvrv, best, sigma, cost, parent, verbose)
            end
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("s: " * string(s), true, false)
        logging("sigma: " * string(sigma), true, false)
        logging("cost: " * string(cost), true, false)
        logging("parent: " * string(parent), true, false)
        logging("====================================================", true, false)
    end
    return sigma, cost, parent
end

# Counting number of optimal walks
# function target_cost(e, c)
#     return e[3]+e[4] - c 
# end

function target_cost(e, c)
    return c
end

function min_cost()
    return 0
end

function max_cost()
    return typemax(Int64)
end

function sharp_values(n, earr, s, sigma, cost, verbose::Bool)
    optimal_value_v = fill(max_cost(), n)
    optimal_value_v[s] = min_cost()
    for ei in 1:lastindex(earr)
        if (!ismissing(cost[ei]))
            c = target_cost(earr[ei], cost[ei])
            if (optimal_value_v[earr[ei][2]] > c)
                optimal_value_v[earr[ei][2]] = c
            end
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("optimal_value_v: " * string(optimal_value_v), true, false)
        logging("====================================================", true, false)
    end
    sharp_e = fill(0, length(earr))
    for ei in 1:lastindex(earr)
        if (!ismissing(cost[ei]))
            c = target_cost(earr[ei], cost[ei])
            if (optimal_value_v[earr[ei][2]] > 0 && c == optimal_value_v[earr[ei][2]])
                sharp_e[ei] = sigma[ei]
            end
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("sharp_e: " * string(sharp_e), true, false)
        logging("====================================================", true, false)
    end
    return sharp_e
end

function sigma_node(n, earr, s, sharp_e, verbose)
    sigma_v = fill(0, n)
    sigma_v[s] = 1
    for ei in 1:lastindex(earr)
        v = earr[ei][2]
        sigma_v[v] = sigma_v[v] + sharp_e[ei]
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("sigma_v: " * string(sigma_v), true, false)
        logging("====================================================", true, false)
    end
    return sigma_v
end

function init_b(earr, sharp, sigma_v, verbose)
    b = zeros(length(earr))
    for ei in 1:lastindex(earr)
        v = earr[ei][2]
        if (sigma_v[v] > 0)
            b[ei] = sharp[ei] / sigma_v[v]
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("b (initial values): " * string(b), true, false)
        logging("====================================================", true, false)
    end
    return b
end

function temporal_walk_betweenness_s(n, alpha, beta, earr, edepv, edepv_index, s, verbose)
    # start_time = time()
    sigma_e, cost_e, parent_e = optimal_walks_counter(n, alpha, beta, earr, edepv, edepv_index, s, false)
    # println(s, " (forward): ", time() - start_time)
    sharp_e = sharp_values(n, earr, s, sigma_e, cost_e, verbose)
    sigma_v = sigma_node(n, earr, s, sharp_e, verbose)
    # start_time = time()
    b_e = init_b(earr, sharp_e, sigma_v, false)
    for ei in lastindex(earr):-1:1
        if (b_e[ei] > 0 && length(parent_e[ei]) > 0)
            for p in parent_e[ei]
                b_e[p] = b_e[p] + b_e[ei] * sigma_e[p] / sigma_e[ei]
            end
        end
    end
    # println(s, " (backward): ", time() - start_time)
    if (verbose)
        logging("====================================================", true, false)
        logging("b_e: " * string(b_e), true, false)
        logging("====================================================", true, false)
    end
    return b_e, sigma_v
end

function temporal_walk_betweenness(fn::String, sep, verbose)
    n, alpha, beta, earr, edep = read_patg(fn, sep, false)
    edepv, edepv_index = distribute_edep(n, earr, edep, verbose)
    b_e = zeros(length(earr))
    b_v = zeros(n)
    for s in 1:n
        b_e_s, sigma_v_s = temporal_walk_betweenness_s(n, alpha, beta, earr, edepv, edepv_index, s, false)
        b_e = b_e .+ b_e_s
        for v in 1:n
            if (s != v && sigma_v_s[v] > 0)
                b_v[v] = b_v[v] - 1
            end
        end
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("b_e: " * string(b_e), true, false)
        logging("====================================================", true, false)
    end
    for ei in 1:lastindex(earr)
        v = earr[ei][2]
        b_v[v] = b_v[v] + b_e[ei]
    end
    if (verbose)
        logging("====================================================", true, false)
        logging("b_v: " * string(b_v), true, false)
        logging("====================================================", true, false)
    end
    return b_v
end
