struct TG
    num_nodes::Int64
    temporal_edges::Array{Tuple{Int64,Int64,Int64}}
    file_id::Array{String}
    file_time::Array{Int64}
end

function read_tg(file_name::String, sep::String)
    current_node_id::Int64 = 1
    file_id::Vector{String} = []
    file_id_to_graph_id::Dict{String,Int64} = Dict{String,Int64}()
    current_time::Int64 = 1
    file_time::Vector{Int64} = []
    file_time_to_graph_time::Dict{Int64,Int64} = Dict{Int64,Int64}()
    temporal_edges::Array{Tuple{Int64,Int64,Int64}} = []
    t::Int64 = 0
    f::IOStream = open(file_name, "r")
    first = true
    for line in eachline(f)
        if (first)
            first = false
        else
            split_line::Vector{String} = split(line, sep)
            t = parse(Int64, split_line[3])
            if (!haskey(file_id_to_graph_id, split_line[1]))
                file_id_to_graph_id[split_line[1]] = current_node_id
                push!(file_id, split_line[1])
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_id_to_graph_id, split_line[2]))
                file_id_to_graph_id[split_line[2]] = current_node_id
                push!(file_id, split_line[2])
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_time_to_graph_time, t))
                file_time_to_graph_time[t] = current_time
                push!(file_time, t)
                current_time = current_time + 1
            end
        end
    end
    close(f)
    sort!(file_time)
    for t in 1:lastindex(file_time)
        file_time_to_graph_time[file_time[t]] = t
    end
    f = open(file_name, "r")
    first = true
    for line in eachline(f)
        if (first)
            first = false
        else
            split_line::Vector{String} = split(line, sep)
            t = parse(Int64, split_line[3])
            push!(temporal_edges, (file_id_to_graph_id[split_line[1]], file_id_to_graph_id[split_line[2]], file_time_to_graph_time[t]))
        end
    end
    sort!(temporal_edges, by=te -> te[3])
    return TG(length(file_id_to_graph_id), temporal_edges, file_id, file_time)
end

function print_tg_stats(tg; graph_name="anonymous")
    logging("====================================================", true, false)
    logging("Temporal network: " * graph_name, true, false)
    logging("====================================================", true, false)
    logging("Number of nodes " * string(tg.num_nodes), true, false)
    logging("Number temporal of edges " * string(length(tg.temporal_edges)), true, false)
    logging("Number of unique time stamps " * string(length(tg.file_time)), true, false)
    logging("====================================================", true, false)
end

function temporal_adjacency_list(tg)
    tal::Array{Array{Tuple{Int64,Int64}}} = Array{Array{Tuple{Int64,Int64}}}(undef, tg.num_nodes)
    for u in 1:tg.num_nodes
        tal[u] = Tuple{Int64,Int64,Int64}[]
    end
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        push!(tal[te[1]], (te[2], te[3]))
    end
    return tal
end

function temporal_node_index(tg)
    d::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    current_index = 1
    for edge in tg.temporal_edges
        if (get(d, (edge[2], edge[3]), 0) == 0)
            d[(edge[2], edge[3])] = current_index
            current_index = current_index + 1
        end
    end
    for s in 1:tg.num_nodes
        d[(s, 0)] = current_index
        current_index = current_index + 1
    end
    return d
end

function next_temporal_neighbours(tn_index::Dict{Tuple{Int64,Int64},Int64}, tal::Array{Array{Tuple{Int64,Int64}}})
    ntn_pos::Array{Int64} = zeros(length(tn_index))
    for tn in keys(tn_index)
        u = tn[1]
        t = tn[2]
        tni = tn_index[tn]
        left::Int64 = 1
        right::Int64 = length(tal[u])
        pos::Int64 = length(tal[u]) + 1
        mid::Int64 = -1
        while (left <= right)
            mid = (left + right) ÷ 2
            if (tal[u][mid][2] <= t)
                left = mid + 1
            else
                pos = mid
                right = mid - 1
            end
        end
        ntn_pos[tni] = pos
    end
    return ntn_pos
end

function bmnr_shortest(fn, sep, verbose_step)
    tg = read_tg(fn, sep)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    ntn_pos::Array{Int64} = next_temporal_neighbours(tn_index, tal)
    nn = tg.num_nodes
    ntn = length(keys(tn_index))
    # Begin data structures for BFS
    sigma::Array{UInt128} = Array{UInt128}(undef, nn)
    dist::Array{Int64} = Array{Int64}(undef, nn)
    sigma_t::Array{UInt128} = Array{UInt128}(undef, ntn)
    dist_t::Array{Int64} = Array{Int64}(undef, ntn)
    delta_sh::Array{Float64} = Array{Float64}(undef, ntn)
    predecessors::Array{Set{Tuple{Int64,Int64}}} = Array{Set{Tuple{Int64,Int64}}}(undef, ntn)
    queue::Queue{Tuple{Int64,Int64}} = Queue{Tuple{Int64,Int64}}()
    stack::Stack{Tuple{Int64,Int64}} = Stack{Tuple{Int64,Int64}}()
    # End data structures for BFS
    temporal_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    v::Int64 = -1
    t_v::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    processed_so_far::Int64 = 0
    start_time = time()
    for s in 1:tg.num_nodes
        for u in 1:tg.num_nodes
            dist[u] = -1
            sigma[u] = 0
        end
        for tn in 1:(length(tn_index))
            sigma_t[tn] = 0
            delta_sh[tn] = 0
            dist_t[tn] = -1
            predecessors[tn] = Set{Tuple{Int64,Int64}}()
        end
        tni = tn_index[(s, 0)]
        sigma[s] = 1
        sigma_t[tni] = 1
        dist[s] = 0
        dist_t[tni] = 0
        enqueue!(queue, (s, 0))
        iter = 0
        while length(queue) != 0
            iter += 1
            temporal_node = dequeue!(queue)
            u = temporal_node[1]
            t = temporal_node[2]
            tni = tn_index[(u, t)]
            for neig in tal[u][ntn_pos[tni]:end]
                w = neig[1]
                t_w = neig[2]
                tni_w = tn_index[(w, t_w)]
                if dist_t[tni_w] == -1
                    dist_t[tni_w] = dist_t[tni] + 1
                    if dist[w] == -1
                        dist[w] = dist_t[tni_w]
                    end
                    enqueue!(queue, neig)
                    push!(stack, neig)
                end
                if dist_t[tni_w] == dist_t[tni] + 1
                    if (sigma_t[tni] > typemax(UInt128) - sigma_t[tni_w])
                        logging("Overflow occurred with source " * s)
                        return [], 0.0
                    end
                    sigma_t[tni_w] += sigma_t[tni]
                    push!(predecessors[tni_w], temporal_node)
                    if dist_t[tni_w] == dist[w]
                        if (sigma_t[tni] > typemax(UInt128) - sigma[w])
                            logging("Overflow occurred with source " * s)
                            return [], 0.0
                        end
                        sigma[w] += sigma_t[tni]
                    end
                end
            end
        end
        temporal_betweenness_centrality[s] -= (count(x -> x >= 0, dist) - 1)
        while length(stack) != 0
            temporal_node = pop!(stack)
            w = temporal_node[1]
            t_w = temporal_node[2]
            tni_w = tn_index[(w, t_w)]
            if dist_t[tni_w] == dist[w]
                delta_sh[tni_w] += sigma_t[tni_w] / sigma[w]
            end
            for pred in predecessors[tni_w]
                v = pred[1]
                t_v = pred[2]
                tni_v = tn_index[(v, t_v)]
                delta_sh[tni_v] += (sigma_t[tni_v] / sigma_t[tni_w]) * delta_sh[tni_w]
                temporal_betweenness_centrality[v] += (sigma_t[tni_v] / sigma_t[tni_w]) * delta_sh[tni_w]
            end
        end
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            logging("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_betweenness_centrality, finish_total
end

function bmnr_shortest_foremost(fn, sep, verbose_step)
    tg = read_tg(fn, sep)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    ntn_pos::Array{Int64} = next_temporal_neighbours(tn_index, tal)
    nn = tg.num_nodes
    ntn = length(keys(tn_index))
    # Begin data structures for BFS
    sigma::Array{UInt128} = Array{UInt128}(undef, nn)
    dist::Array{Int64} = Array{Int64}(undef, nn)
    sigma_t::Array{UInt128} = Array{UInt128}(undef, ntn)
    dist_t::Array{Int64} = Array{Int64}(undef, ntn)
    predecessors::Array{Set{Tuple{Int64,Int64}}} = Array{Set{Tuple{Int64,Int64}}}(undef, ntn)
    queue::Queue{Tuple{Int64,Int64}} = Queue{Tuple{Int64,Int64}}()
    stack::Stack{Tuple{Int64,Int64}} = Stack{Tuple{Int64,Int64}}()
    t_min::Array{Int64} = Array{Int64}(undef, nn)
    delta_fm::Array{Float64} = Array{Float64}(undef, ntn)
    # End data structures for BFS
    temporal_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    v::Int64 = -1
    t_v::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    processed_so_far::Int64 = 0
    start_time = time()
    for s in 1:tg.num_nodes
        for u in 1:tg.num_nodes
            dist[u] = -1
            sigma[u] = 0
            t_min[u] = -1
        end
        for tn in 1:(length(tn_index))
            sigma_t[tn] = 0
            dist_t[tn] = -1
            predecessors[tn] = Set{Tuple{Int64,Int64}}()
            delta_fm[tn] = 0
        end
        tni = tn_index[(s, 0)]
        sigma[s] = 1
        t_min[s] = 0
        sigma_t[tni] = 1
        dist[s] = 0
        dist_t[tni] = 0
        enqueue!(queue, (s, 0))
        iter = 0
        while length(queue) != 0
            iter += 1
            temporal_node = dequeue!(queue)
            u = temporal_node[1]
            t = temporal_node[2]
            tni = tn_index[(u, t)]
            for neig in tal[u][ntn_pos[tni]:end]
                w = neig[1]
                t_w = neig[2]
                tni_w = tn_index[(w, t_w)]
                if dist_t[tni_w] == -1
                    dist_t[tni_w] = dist_t[tni] + 1
                    if dist[w] == -1
                        dist[w] = dist_t[tni_w]
                    end
                    enqueue!(queue, neig)
                    push!(stack, neig)
                end
                if dist_t[tni_w] == dist_t[tni] + 1
                    if (sigma_t[tni] > typemax(UInt128) - sigma_t[tni_w])
                        logging("Overflow occurred with source " * s)
                        return [], 0.0
                    end
                    sigma_t[tni_w] += sigma_t[tni]
                    push!(predecessors[tni_w], temporal_node)
                    if dist_t[tni_w] == dist[w]
                        if (sigma_t[tni] > typemax(UInt128) - sigma[w])
                            logging("Overflow occurred with source " * s)
                            return [], 0.0
                        end
                        sigma[w] += sigma_t[tni]
                    end
                end
                if (t_min[w] == -1 || t_w < t_min[w])
                    t_min[w] = t_w
                end
            end
        end
        temporal_betweenness_centrality[s] -= (count(x -> x >= 0, dist) - 1)
        while length(stack) != 0
            temporal_node = pop!(stack)
            w = temporal_node[1]
            t_w = temporal_node[2]
            tni_w = tn_index[(w, t_w)]
            if (t_min[w] == t_w)
                delta_fm[tni_w] += 1
            end
            for pred in predecessors[tni_w]
                v = pred[1]
                t_v = pred[2]
                tni_v = tn_index[(v, t_v)]
                delta_fm[tni_v] += (sigma_t[tni_v] / sigma_t[tni_w]) * delta_fm[tni_w]
                temporal_betweenness_centrality[v] += (sigma_t[tni_v] / sigma_t[tni_w]) * delta_fm[tni_w]
            end
        end
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            logging("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_betweenness_centrality, finish_total
end
