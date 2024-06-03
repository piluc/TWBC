@inline function finalize_costs_nr_shortest_foremost(u, ie, patg, l, c, σ, C, Σ)
    if (c[u] < typemax(Int64))
        for j in max(l[u], 1):ie
            f = patg.edepv[u][j]
            C[f] = c[u] + 1
            Σ[f] = σ[u]
        end
    end
    l[u] = ie + 1
end

function nrsfo(fn, sep, verbose_step)
    patg = read_patg(fn, sep, α=0, β=typemax(Int64))
    n::Int64 = patg.n
    M::Int64 = length(patg.earr)
    b_v::Array{Float64} = zeros(Float64, n)
    # Begin data structures for BCV algorithm
    l::Array{Int64} = [1 for _ = 1:n]
    c::Array{Int64} = [typemax(Int64) for _ = 1:n]
    σ::Array{UInt128} = [0 for _ = 1:n]
    σstar::Array{UInt128} = [0 for _ = 1:n]
    # cstar::Array{Array{Int64}} = [[typemax(Int64), typemax(Int64)] for _ = 1:n]
    cstar::Matrix{Int64} = fill(typemax(Int64), n, 2)
    δ::Array{Float64} = [0 for _ = 1:n]
    L::Array{Int64} = [0 for _ = 1:M]
    C::Array{Int64} = [typemax(Int64) for _ = 1:M]
    Σ::Array{UInt128} = [0 for _ = 1:M]
    Σstar::Array{UInt128} = [0 for _ = 1:M]
    b::Array{Float64} = [0 for _ = 1:M]
    # End data structures for BCV algorithm
    processed_so_far::Int64 = 0
    start_time = time()
    for s in 1:n
        # Begin forward phase
        for ei in 1:M
            e = patg.earr[ei]
            u, v, τ, λ = e
            ie = patg.edepv_index[ei]
            if (ie >= l[u])
                finalize_costs_nr_shortest_foremost(u, ie, patg, l, c, σ, C, Σ)
            end
            if (u == s)
                C[ei] = 1
                Σ[ei] = 1
            end
            if (C[ei] != typemax(Int64))
                if (C[ei] <= c[v])
                    a = l[v]
                    while (a <= length(patg.edepv[v]) && patg.earr[patg.edepv[v][a]][3] < τ + λ)
                        a = a + 1
                    end
                    finalize_costs_nr_shortest_foremost(v, a - 1, patg, l, c, σ, C, Σ)
                    if (C[ei] < c[v])
                        c[v] = C[ei]
                        σ[v] = 0
                    end
                    if (Σ[ei] > typemax(UInt128) - σ[v])
                        logging("Overflow have occurred with source " * string(s), true, false)
                        logging("Trying to add " * string(Σ[ei]) * " to " * string(σ[v]), true, false)
                        return []
                    end
                    σ[v] += Σ[ei]
                    L[ei] = a
                end
            end
        end
        # End forward phase
        # Begin pre-backward phase
        for ei in 1:M
            if (C[ei] < typemax(Int64))
                _, v, τ, λ = patg.earr[ei]
                if ((τ + λ < cstar[v, 1]) || (τ + λ == cstar[v, 1] && C[ei] < cstar[v, 2]))
                    cstar[v, 1] = τ + λ
                    cstar[v, 2] = C[ei]
                end
            end
        end
        # return cstar
        for ei in 1:M
            if (C[ei] < typemax(Int64))
                _, v, τ, λ = patg.earr[ei]
                if (τ + λ == cstar[v, 1] && C[ei] == cstar[v, 2])
                    Σstar[ei] = Σ[ei]
                    σstar[v] += Σ[ei]
                end
            end
        end
        # End pre-backward phase
        # Begin backward phase
        for ei in M:-1:1
            _, v, τ, λ = patg.earr[ei]
            if (v != s && L[ei] > 0)
                if (τ + λ > cstar[v, 1] || cstar[v, 2] < C[ei])
                    δ[v] = 0.0
                    cstar[v, 1] = τ + λ
                    cstar[v, 2] = C[ei]
                end
                # if (c[v] < C[ei])
                #     δ[v] = 0.0
                #     c[v] = C[ei]
                # end
                for i in L[ei]:(l[v]-1)
                    fi = patg.edepv[v][i]
                    δ[v] += b[fi] / Σ[fi]
                end
                l[v] = L[ei]
                b[ei] = Σ[ei] * δ[v]
                if (Σstar[ei] > 0)
                    b[ei] += Σstar[ei] / σstar[v]
                end
            end
            # println(string(patg.earr[ei]), " \\idep{", L[ei], "} & \$", C[ei], "\$ & \$", Σ[ei], "\$ & \$", Σstar[ei], "\$ & \$", b[ei], "\$ & \$", string(patg.edepv[v]), "\$ & \\idep{", l[v], "} & \$", c[v], "\$ & \$", σstar[v], "\$ & \$", δ[v], "\$")
        end
        # End backward phase
        # Accumulate edge betweenness in nodes
        b_v_s::Array{Float64} = zeros(Float64, n)
        for ei in 1:M
            v = patg.earr[ei][2]
            if (s != v)
                b_v_s[v] += b[ei]
            end
        end
        for v in 1:n
            if (s != v && σstar[v] > 0)
                b_v_s[v] -= 1
            end
            b_v[v] += b_v_s[v]
        end
        # Re-initialise data structures for next source
        l = [1 for _ = 1:n]
        c = [typemax(Int64) for _ = 1:n]
        σ = [0 for _ = 1:n]
        σstar = [0 for _ = 1:n]
        # cstar = [[typemax(Int64), typemax(Int64)] for _ = 1:n]
        cstar = fill(typemax(Int64), n, 2)
        δ = [0 for _ = 1:n]
        L = [0 for _ = 1:M]
        C = [typemax(Int64) for _ = 1:M]
        Σ = [0 for _ = 1:M]
        Σstar = [0 for _ = 1:M]
        b = [0 for _ = 1:M]
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            logging("TSB. Processed " * string(processed_so_far) * "/" * string(n) * " nodes in " * finish_partial * " seconds")
        end
    end
    return b_v, time() - start_time
end
