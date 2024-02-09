LargeInt = UInt128;
LargeReal = Float64;
LargeExact = false;
# LargeInt = BigInt; LargeReal = Float64; LargeExact = false;
# LargeInt = BigInt; LargeReal = BigFloat; LargeExact = false;
# LargeInt = BigInt; LargeReal = BigRational; LargeExact = true;

mutable struct I_ITEM_DH
    l::Int64
    r::Int64
    c::DH_COST
    P::Vector{Int64}
    η::LargeInt
    function I_ITEM_DH(_l, _r, _c, _P, _η)
        return new(_l, _r, _c, _P, _η)
    end
end

struct FWDBWD_5
    I::Array{Vector{I_ITEM_DH}}
    l::Array{Int64}
    l_refresh::Array{Int64}
    r::Array{Int64}
    γ::Array{DH_COST}
    C::Array{Union{DH_COST,Missing}}
    σe::Array{LargeInt}
    σe_star::Array{LargeInt}
    σv_star::Array{LargeInt}
    c_star::Array{Union{TARGET_COST,Missing}}
    be::Array{LargeReal}
    bv::Array{LargeReal}
    Left::Array{Int64}
    Right::Array{Int64}
    function FWDBWD_5(tg, ct::CS)
        n = tg.n
        M = length(tg.earr)
        _I::Array{Vector{I_ITEM_DH}} = [Vector{I_ITEM_DH}() for _ = 1:n]
        _l::Array{Int64} = fill(1, n)
        _l_refresh::Array{Int64} = [length(tg.edepv[v]) + 1 for v = 1:n]
        _r::Array{Int64} = fill(0, n)
        _γ::Array{DH_COST} = [cs_edge_cost(ct, tg.earr[ei]) for ei = 1:M]
        _C = Array{Union{DH_COST,Missing}}(missing, M)
        _σe::Array{LargeInt} = fill(0, M)
        _σe_star::Array{LargeInt} = fill(0, M)
        _σv_star::Array{LargeInt} = fill(0, n)
        _c_star = Array{Union{TARGET_COST,Missing}}(missing, n)
        _be::Array{LargeReal} = fill(0, M)
        _bv::Array{LargeReal} = fill(0, n)
        _Left::Array{Int64} = fill(M + 1, M)
        _Right::Array{Int64} = fill(0, M)
        return new(_I, _l, _l_refresh, _r, _γ, _C, _σe, _σe_star, _σv_star, _c_star, _be, _bv, _Left, _Right)
    end
end

@inline function finalize_cost(v, j, tg, ct, fb5)
    Iv = fb5.I[v]
    while length(Iv) > 0 && first(Iv).l <= j
        first_I::I_ITEM_DH = first(Iv)
        l = first_I.l
        P = first_I.P
        η = first_I.η
        while length(P) > 0 && fb5.Right[first(P)] <= j
            ei = popfirst!(P)
            for f in tg.edepv[v][l:fb5.Right[ei]]
                fb5.C[f] = cs_oplus(ct, first_I.c, fb5.γ[f])
                fb5.σe[f] = η
            end
            η -= fb5.σe[ei]
            l = fb5.Right[ei] + 1
        end
        for f in tg.edepv[v][l:min(j, first_I.r)]
            fb5.C[f] = cs_oplus(ct, first_I.c, fb5.γ[f])
            fb5.σe[f] = η
        end
        if j >= first_I.r
            popfirst!(Iv)
        else
            first_I.l = j + 1
            first_I.η = η
            break
        end
    end
    fb5.l[v] = j + 1
end

function forward_phase(tg, s, ct, fb5::FWDBWD_5)
    n::Int64 = tg.n
    M::Int64 = length(tg.earr)
    for ei in 1:M
        e = tg.earr[ei]
        u = e[1]
        v = e[2]
        arr = e[3] + e[4]
        i = tg.edepv_index[ei]
        if i >= fb5.l[u]
            finalize_cost(u, i, tg, ct, fb5)
        end
        if u == s || !ismissing(fb5.C[ei])
            if u == s
                if ismissing(fb5.C[ei]) || cs_lessthan(ct, fb5.γ[ei], fb5.C[ei])
                    fb5.C[ei] = fb5.γ[ei]
                    fb5.σe[ei] = 1
                elseif cs_equal(ct, fb5.γ[ei], fb5.C[ei])
                    if ((LargeInt != BigInt) && fb5.σe[ei] == typemax(LargeInt))
                        logging("Overflow have occurred with source " * string(s), true, false)
                        logging("Trying to add 1 to " * string(fb5.σe[ei]), true, false)
                        return []
                    end
                    fb5.σe[ei] += 1
                end
            end
            arr_α = arr + tg.α[v]
            β = tg.β[v]
            l = fb5.l[v]
            while (l <= length(tg.edepv[v]) && tg.edepv_dep[v][l] < arr_α)
                l += 1
            end
            r = fb5.r[v]
            while (r + 1 <= length(tg.edepv[v]) && tg.edepv_dep[v][r+1] - arr <= β)
                r += 1
            end
            finalize_cost(v, l - 1, tg, ct, fb5)
            l_c = max(l, fb5.r[v] + 1)
            Iv = fb5.I[v]
            while length(Iv) > 0 && cs_lessthan(ct, fb5.C[ei], last(Iv).c)
                for ej in last(Iv).P
                    fb5.Right[ej] = l - 1
                    l_c = last(Iv).l
                end
                pop!(Iv)
            end
            if length(Iv) > 0 && cs_equal(ct, fb5.C[ei], last(Iv).c)
                last_I = last(Iv)
                last_I.r = r
                if ((LargeInt != BigInt) && fb5.σe[ei] > typemax(LargeInt) - last_I.η)
                    logging("Overflow have occurred with source " * string(s), true, false)
                    logging("Trying to add " * string(fb5.σe[ei]) * " to " * string(last_I.η), true, false)
                    return []
                end
                last_I.η += fb5.σe[ei]
                push!(last_I.P, ei)
                fb5.Left[ei] = last_I.l
            elseif l_c <= r
                fb5.Left[ei] = l_c
                push!(fb5.I[v], I_ITEM_DH(l_c, r, fb5.C[ei], [ei], fb5.σe[ei]))
            end
            fb5.Right[ei] = r
            fb5.r[v] = r
        end
    end
end

function backward_phase(tg, s, ct, fb5::FWDBWD_5)
    n::Int64 = tg.n
    M::Int64 = length(tg.earr)
    fb5.c_star[s] = cs_min_target_cost(ct)
    for ei = 1:M
        if !ismissing(fb5.C[ei])
            v = tg.earr[ei][2]
            c = cs_target_cost(ct, tg.earr[ei], fb5.C[ei])
            if ismissing(fb5.c_star[v]) || cs_target_lessthan(ct, c, fb5.c_star[v])
                fb5.c_star[v] = c
            end
        end
    end
    fb5.σv_star[s] = 1
    for ei = 1:M
        if !ismissing(fb5.C[ei])
            v = tg.earr[ei][2]
            c = cs_target_cost(ct, tg.earr[ei], fb5.C[ei])
            if cs_target_equal(ct, c, fb5.c_star[v])
                fb5.σe_star[ei] = fb5.σe[ei]
                fb5.σv_star[v] += fb5.σe[ei]
            end
        end
    end
    for ei = M:-1:1
        v = tg.earr[ei][2]
        if fb5.Left[ei] <= fb5.Right[ei]
            lv = fb5.l[v]
            for f in tg.edepv[v][max(lv, fb5.Right[ei] + 1):fb5.r[v]]
                if v != s || cs_equal(ct, fb5.C[f], cs_oplus(ct, fb5.C[ei], fb5.γ[f]))
                    fb5.bv[v] -= fb5.be[f] / fb5.σe[f]
                end
            end
            fb5.r[v] = fb5.Right[ei]
            rv = fb5.r[v]
            if (!LargeExact) && (rv < lv || rv < fb5.l_refresh[v])
                fb5.bv[v] = 0
                fb5.l[v] = rv + 1
                lv = fb5.l[v]
                fb5.l_refresh[v] = fb5.Left[ei]
            end
            for f in tg.edepv[v][fb5.Left[ei]:min(fb5.Right[ei], lv - 1)]
                if v != s || cs_equal(ct, fb5.C[f], cs_oplus(ct, fb5.C[ei], fb5.γ[f]))
                    fb5.bv[v] += fb5.be[f] / fb5.σe[f]
                end
            end
            fb5.l[v] = fb5.Left[ei]
            fb5.be[ei] = fb5.σe[ei] * fb5.bv[v]
        end
        if fb5.σe_star[ei] > 0
            fb5.be[ei] += LargeReal(fb5.σe_star[ei]) / LargeReal(fb5.σv_star[v])
        end
    end
end

function s_temporal_betweenness_5(tg, s, ct, fb5::FWDBWD_5)
    forward_phase(tg, s, ct, fb5)
    backward_phase(tg, s, ct, fb5)
end

function temporal_betweenness_5(fn::String, sep, ct::CS, v_step::Int64; _α=0, _β=typemax(Int64))
    tg = read_patg(fn, sep, α=_α, β=_β)
    n::Int64 = tg.n
    M::Int64 = length(tg.earr)
    b_v::Array{LargeReal} = zeros(LargeReal, n)
    processed_so_far::Int64 = 0
    start_time = time()
    for s in 1:n
        fb5 = FWDBWD_5(tg, ct)
        s_temporal_betweenness_5(tg, s, ct, fb5)
        b_v_s::Array{LargeReal} = zeros(LargeReal, n)
        for ei in 1:M
            v = tg.earr[ei][2]
            if s != v
                b_v_s[v] += fb5.be[ei]
            end
        end
        for v in 1:n
            if (s != v && fb5.σv_star[v] > 0)
                b_v_s[v] -= LargeReal(1)
            end
            b_v[v] += b_v_s[v]
        end
        processed_so_far = processed_so_far + 1
        if (v_step > 0 && processed_so_far % v_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            logging("Processed " * string(processed_so_far) * "/" * string(n) * " nodes in " * finish_partial * " seconds")
        end
    end
    return b_v, time() - start_time
end

function tfab(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TFaB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tfob(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TFoB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tlab(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TLaB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tsb(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TSB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tsfab(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TSFaB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tsfob(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TSFoB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end

function tslab(fn::String, sep, v_step; α=0, β=typemax(Int64))
    cs = TSLaB_CT()
    return temporal_betweenness_5(fn, sep, cs, v_step, _α=α, _β=β)
end