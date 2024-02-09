mutable struct I_ITEM
    l::Int64
    r::Int64
    c::Int64
    P::Vector{Int64}
    η::UInt128
    function I_ITEM(_l, _r, _c, _P, _η)
        return new(_l, _r, _c, _P, _η)
    end
end

@inline function finalize_cost_shortest_foremost(v, j, tg, I, l, R, C, Σ)
    while (length(I[v]) > 0 && first(I[v]).l <= j)
        Q = first(I[v])
        ql = Q.l
        qη = Q.η
        while (length(Q.P) > 0 && R[first(Q.P)] <= j)
            ei = popfirst!(Q.P)
            for fi in tg.edepv[v][ql:R[ei]]
                C[fi] = Q.c + 1
                Σ[fi] = qη
            end
            qη = qη - Σ[ei]
            ql = R[ei] + 1
            @assert ql <= Q.r + 1 "buggy Right?"
        end
        for fi in tg.edepv[v][ql:min(j, Q.r)]
            C[fi] = Q.c + 1
            Σ[fi] = qη
        end
        if (j >= Q.r)
            popfirst!(I[v])
        else
            Q.l = j + 1
            Q.η = qη
            break
        end
    end
    l[v] = j + 1
end

function algorithm3(fn, sep, verbose_step; _α=0, _β=typemax(Int64))
    patg = read_patg(fn, sep, α=_α, β=_β)
    n::Int64 = patg.n
    M::Int64 = length(patg.earr)
    # Begin data structures for BCV algorithm
    l::Array{Int64} = ones(n)
    r::Array{Int64} = zeros(n)
    I::Array{Vector{I_ITEM}} = [Vector{I_ITEM}() for _ = 1:n]
    L::Array{Int64} = [length(patg.earr) + 1 for _ = 1:M]
    R::Array{Int64} = [0 for _ = 1:M]
    C::Array{Int64} = [typemax(Int64) for _ = 1:M]
    Σ::Array{UInt128} = [0 for _ = 1:M]
    cstar::Matrix{Int64} = fill(typemax(Int64), n, 2)
    σstar::Array{UInt128} = zeros(UInt128, n)
    δ::Array{Float64} = [0 for _ = 1:n]
    Σstar::Array{UInt128} = zeros(UInt128, M)
    b_e::Array{Float64} = [0 for _ = 1:M]
    l_refresh::Array{Int64} = [length(patg.edepv[v]) + 1 for v = 1:n]
    # End data structures for BCV algorithm
    b_v::Array{Float64} = zeros(Float64, n)
    processed_so_far::Int64 = 0
    start_time = time()
    for s in 1:n
        # Begin forward phase
        for ei in 1:M
            u, v, τ, λ = patg.earr[ei]
            i = patg.edepv_index[ei]
            if (i >= l[u])
                finalize_cost_shortest_foremost(u, i, patg, I, l, R, C, Σ)
            end
            if (u == s)
                C[ei] = 1
                Σ[ei] = 1
            end
            if (C[ei] < typemax(Int64))
                a = l[v]
                while (a <= length(patg.edepv[v]) && patg.earr[patg.edepv[v][a]][3] < τ + λ)
                    a = a + 1
                end
                b = r[v]
                while (b < length(patg.edepv[v]) && patg.earr[patg.edepv[v][b+1]][3] - (τ + λ) <= _β)
                    b = b + 1
                end
                finalize_cost_shortest_foremost(v, a - 1, patg, I, l, R, C, Σ)
                lc = max(a, r[v] + 1)
                while (length(I[v]) > 0 && C[ei] < last(I[v]).c)
                    Q = pop!(I[v])
                    lc = Q.l
                    for fi in Q.P
                        R[fi] = a - 1
                    end
                end
                if (length(I[v]) > 0 && C[ei] == last(I[v]).c)
                    Q = last(I[v])
                    Q.r = b
                    if (Σ[ei] > typemax(UInt128) - Q.η)
                        logging("Overflow have occurred with source " * string(s), true, false)
                        logging("Trying to add " * string(Σ[ei]) * " to " * string(Q.η), true, false)
                        return []
                    end
                    Q.η += Σ[ei]
                    push!(Q.P, ei)
                    L[ei] = Q.l
                elseif (lc <= b)
                    L[ei] = lc
                    push!(I[v], I_ITEM(lc, b, C[ei], [ei], Σ[ei]))
                end
                R[ei] = b
                r[v] = b
            end
        end
        # End forward phase
        # Begin pre-backward phase
        cstar[s, 1] = 0
        cstar[s, 2] = 0
        σstar[s] = 1
        for ei in 1:M
            _, v, τ, λ = patg.earr[ei]
            if (C[ei] < typemax(Int64) && ((τ + λ < cstar[v, 1]) || (τ + λ == cstar[v, 1] && C[ei] < cstar[v, 2])))
                cstar[v, 1] = τ + λ
                cstar[v, 2] = C[ei]
            end
        end
        for ei in 1:M
            _, v, τ, λ = patg.earr[ei]
            if (τ + λ == cstar[v, 1] && C[ei] == cstar[v, 2])
                Σstar[ei] = Σ[ei]
                if (σstar[v] > typemax(UInt128) - Σ[ei])
                    logging("Overflow have occurred with source " * string(s), true, false)
                    logging("Trying to add " * string(σstar[v]) * " to " * string(Σ[ei]), true, false)
                    return []
                end
                σstar[v] += Σ[ei]
            end
        end
        # End pre-backward phase
        # Begin backward phase
        for ei in M:-1:1
            _, v, _, _ = patg.earr[ei]
            if (L[ei] <= R[ei])
                for fi in patg.edepv[v][max(l[v], R[ei] + 1):r[v]]
                    if (v != s || C[fi] == C[ei] + 1)
                        δ[v] -= b_e[fi] / Σ[fi]
                    end
                end
                r[v] = R[ei]
                if (r[v] < l[v] || r[v] < l_refresh[v])
                    δ[v] = 0
                    l[v] = r[v] + 1
                    l_refresh[v] = L[ei]
                end
                for fi in patg.edepv[v][L[ei]:min(R[ei], l[v] - 1)]
                    if (v != s || C[fi] == C[ei] + 1)
                        δ[v] += b_e[fi] / Σ[fi]
                    end
                end
                l[v] = L[ei]
                b_e[ei] = Σ[ei] * δ[v]
            end
            if (Σstar[ei] > 0)
                b_e[ei] += Σstar[ei] / σstar[v]
            end
        end
        # End backward phase
        # Accumulate edge betweenness in nodes
        b_v_s::Array{Float64} = zeros(Float64, n)
        for ei in 1:M
            v = patg.earr[ei][2]
            if (s != v)
                b_v_s[v] += b_e[ei]
            end
        end
        for v in 1:n
            if (s != v && σstar[v] > 0)
                b_v_s[v] -= 1
            end
            b_v[v] += b_v_s[v]
        end
        # Re-initialise data structures for next source
        l = ones(n)
        r = zeros(n)
        I = [Vector{I_ITEM}() for _ = 1:n]
        L = [length(patg.earr) + 1 for _ = 1:M]
        R = [0 for _ = 1:M]
        C = [typemax(Int64) for _ = 1:M]
        Σ = [0 for _ = 1:M]
        cstar = fill(typemax(Int64), n, 2)
        σstar = zeros(UInt128, n)
        δ = [0 for _ = 1:n]
        Σstar = zeros(UInt128, M)
        b_e = [0 for _ = 1:M]
        l_refresh = [length(patg.edepv[v]) + 1 for v = 1:n]
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            logging("TSB. Processed " * string(processed_so_far) * "/" * string(n) * " nodes in " * finish_partial * " seconds")
        end
    end
    return b_v, time() - start_time
end
