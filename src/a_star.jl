#= A* Search script that contains the navigation functions.
=#
function a_star_search(start::Int64,
    goal::Int64,
    transitions::Dict{Int64, Vector{LVTransition}},
    heuristic::Function)
    #= Perform an A* search over adjacency lists
		
	Args:
		start: Start state index
        goal: Goal state index
        transitions: Dict{Int64 => Vector{LVTransition}} object that contains 
				all transitions information
        heuristic: Heuristic function that takes in args (index1, index2)

	Return: (tuple)
		success: Boolean for whether the search was successful
        total_path: Total path from start to current, in chronological order
        total_action: Total actions taken, in chronological order
        cost: Cost of the path
	=#

    # Set up data structures
    frontier = PriorityQueue{Int64, Float64}()
    came_from = Dict{Int64, Tuple{Int64, Int64}}()
    graph_costs = Dict{Int64, Float64}()
    full_costs =  Dict{Int64, Float64}()

    # Initialize
    frontier[start] = 0.0
    came_from[start] = (-1, -1)
    graph_costs[start] = 0.0
    full_costs[start] = 0.0
    current = (start, 0.0)

    # Attempt to find path to goal
    while !isempty(frontier)
        # Get the current node
        current = peek(frontier).first
        if current == goal
            return _reconstruct_path(current, goal, came_from, graph_costs)
        end
        delete!(frontier, current) 

        # Iterate over all neighbors (if exists)
        if haskey(transitions, current)
            for next_node in transitions[current]
                next_cost = graph_costs[current] + next_node.cost
                # If this is a better path, update
                if !haskey(graph_costs, next_node.f) || next_cost < graph_costs[next_node.f]
                    came_from[next_node.f] = (current, next_node.a)
                    graph_costs[next_node.f] = next_cost
                    full_costs[next_node.f] = next_cost + heuristic(current, next_node.f)
                    if !haskey(frontier, next_node.f)
                        frontier[next_node.f] = full_costs[next_node.f]
                    end
                end
            end
        end
    end

    # Return failure
    return _reconstruct_path(current, goal, came_from, graph_costs)
end

function _reconstruct_path(current::Int64, 
    goal::Int64,
    came_from::Dict{Int64, Tuple{Int64, Int64}},
    graph_costs::Dict{Int64, Float64})
    #= Reconstruct the result of A* search
		
	Args:
		current: Last ended state index
        goal: Goal state index
        came_from: Dictionary that indicates how the key was obtained,
                showing (prev_index, action_index)
        graph_costs: Dictionary that indicates the cost 
                from the start to state index

	Return: (tuple)
		success: Boolean for whether the search was successful
        total_path: Total path from start to current, in chronological order
        total_action: Total actions taken, in chronological order
        cost: Cost of the path
	=#
    
    # Get path and actions taken
    final = current
    total_path = [current]
    total_action = Int64[]

    # Iterate until we hit the end
    while haskey(came_from, current)
        (current, action) = came_from[current]
        # Append everything but the null pre-start node
        if current > 0 && action > 0
            push!(total_path, current)
            push!(total_action, action)
        end
    end
    reverse!(total_path)
    reverse!(total_action)

    return ((final == goal), total_path, total_action, graph_costs[final])
end


function visualize_a_star_results(n_species::Int64, 
    n_temps::Int64, 
    total_path_list::Vector{Vector{Int64}},
    total_action_list::Vector{Vector{Int64}},
    total_cost_list::Vector{Float64},
    total_pair_list::Vector{Vector{Int64}},
    total_success_list::Vector{Bool},
	candidate_states::Vector{Int64},
	assemblage::Dict{Int64,LVState},
    transitions::Dict{Int64, Vector{LVTransition}},
    perturbations::PerturbationParams,
    costs::TransitionCosts)
    #= Reconstruct the result of A* search
		
	Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
		total_path_list: List of paths from start to current, in chronological order
        total_action_list: List of actions taken, in chronological order
        total_cost_list: List of all costs
        total_pair_list: List of all pairs (start, goal)
        total_success_list: List of all successes
		candidate_states: Vector of all candidate states
        assemblage: Assemblage dictionary
        transitions: Dict{Int64 => Vector{LVTransition}} object that contains 
				all transitions information
        perturbations: Perturbation parameters
        costs: Action costs
        
	Return: 
		a_star_results: DataFrame of all the A star experiment a_star_results
        transition_results: DataFrame of transition statistics
	=#
    # Make data frames for this thread
	a_star_results = DataFrame(
		start = String[],
        goal = String[],
        net_add = Union{Missing, Int64}[],
        net_del = Union{Missing, Int64}[],
        net_temp = Union{Missing, Int64}[],
        full_path = String[],
        path = String[],
        action = String[],
        a_star_cost = Union{Missing, Float64}[],
        a_star_ops = Union{Missing, Int64}[],
        a_star_add = Union{Missing, Int64}[],
        a_star_del = Union{Missing, Int64}[],
        a_star_temp = Union{Missing, Int64}[],)
    transition_results = DataFrame(
        node = String[],
        intermediate_visits = Int64[],
        in_degree = Int64[],
        out_degree = Int64[],)

    # Transition statistics
	transition_visits = Dict{Int64, Int64}()
    transition_in_degree = Dict{Int64, Int64}()
    transition_out_degree = Dict{Int64, Int64}()

    # Get all the path information
    for (total_path, total_action, cost, success, pair) in zip(
        total_path_list, total_action_list, total_cost_list, 
        total_success_list, total_pair_list)
        # Get start and goal
        start = assemblage[pair[1]].str
        goal = assemblage[pair[2]].str

        if success
            # Gather the state evolution list
            path_str_list = []
            for i in total_path
                push!(path_str_list, assemblage[i].str)
                transition_visits[i] = get(transition_visits, i, 0) + 1
            end
            
            # Decrement transition_visits for start and goal
            # we only want intermediate visit frequencies
            transition_visits[total_path[1]] -= 1
            transition_visits[total_path[end]] -= 1

            # Get net changes
            (net_add, net_del, net_temp) = net_species_change(
                total_path[1], total_path[end], assemblage)

            # Get A* total actions
            a_star_add = 0
            a_star_del = 0
            a_star_temp = 0

            # Gather the applied actions list
            action_str_list = []
            for j in total_action
                action = action_idx_to_lv_action(n_species, n_temps, j)
                push!(action_str_list, action.str)
                a_star_add += length(findall(x->x>0, action.ds))
                a_star_del += length(findall(x->x<0, action.ds))
            end
            ops = length(action_str_list)

            # Iterate over action string list (always one less than state list)
            full_path = "[" * start * "]"
            path =  "[" * start * "]"
            action = ""
            for i = 1:length(action_str_list)
                # Manually insert temperature change action (haven't found a better way)
                temp_action = ""
                if path_str_list[i][end] != path_str_list[i+1][end]
                    temp_action = "(*" * path_str_list[i+1][end] * ")"
                    a_star_temp += 1
                end
                action_node = "{" * action_str_list[i] * temp_action * "}"
                next_state_node =  "[" * path_str_list[i+1] * "]"

                full_path *= "~~" * action_node * "~~>" * next_state_node
                path *= "~>" * next_state_node
                action *= action_node * "~"
            end
            push!(a_star_results, [
                start, goal, net_add, net_del, net_temp,
                full_path, path, action, cost, ops, 
                a_star_add, a_star_del, a_star_temp])
        else
            push!(a_star_results, [
                start, goal, missing, missing, missing,
                "N/A", "N/A", "N/A", missing, missing, 
                missing, missing, missing])
        end
    end

    # Get in and out degrees
    for (node, node_transitions) in transitions
        for next_node in node_transitions
            transition_in_degree[next_node.f] = get(transition_in_degree, next_node.f, 0) + 1
            transition_out_degree[node] = get(transition_out_degree, node, 0) + 1
        end
    end

    # Make transition statistics
    for candidate in candidate_states
        push!(transition_results, [
            assemblage[candidate].str, 
            get(transition_visits, candidate, 0),
            get(transition_in_degree, candidate, 0),
            get(transition_out_degree, candidate, 0)])
    end

	return a_star_results, transition_results
end

function optimistic_distance(current::Int64, 
    goal::Int64, 
    assemblage::Dict{Int64,LVState},
    costs::TransitionCosts)
    #= Get optimistic distance heuristic between state 1 and 2. 

    Optimistic distance is the heuristic cost required for A*
    in order to transition from state 1 to state 2, assuming 
    adding all the necessary species can take care of all 
    additions and deletions necessary.
        
    Args:
        start: Start state index
        goal: Goal state index
        assemblage: Assemblage dictionary
        costs: Action costs
    
    Return:
        Float object of optimistic distance
    =#
    # Get the LVStates
    current_state = assemblage[current]
    goal_state = assemblage[goal]

    # Get state distance
    state_dist = sum(max.(goal_state.state - current_state.state, 0))

    # Get temperature distance
    temp_dist = Int64(goal_state.temp != current_state.temp)
    
    return state_dist * costs.add_cost + temp_dist * costs.temp_cost
end

function operation_distance(current::Int64, 
    goal::Int64, 
    assemblage::Dict{Int64,LVState},
    costs::TransitionCosts)
    #= Get operation distance between state 1 and 2. 

    Operation distance is defined by the difference between
    species abundance and temperature between state 1 to state 2.
        
    Args:
        start: Start state index
        goal: Goal state index
        assemblage: Assemblage dictionary
        costs: Action costs
    
    Return:
        Int object of operation distance
    =#
    # Get the LVStates
    current_state = assemblage[current]
    goal_state = assemblage[goal]

    # Get state distance of abundances
    state_dist = Int64(sum(abs.(goal_state.state - current_state.state)))
    
    # Get temperature distance
    temp_dist = Int64(goal_state.temp != current_state.temp)
    
    return state_dist + temp_dist
end

function state_distance(current::Int64, 
    goal::Int64, 
    assemblage::Dict{Int64,LVState},
    perturbations::PerturbationParams,
    costs::TransitionCosts)
    #= Get state distance between state 1 and 2. 

    State distance is defined by the minimum cost required in order
    to transition from state 1 to state 2, assuming no natural transition.
        
    Args:
        start: Start state index
        goal: Goal state index
        assemblage: Assemblage dictionary
        perturbations: Perturbation parameters
        costs: Action costs
    
    Return:
        Float object of state distance
    =#
    # Get the LVStates
    current_state = assemblage[current]
    goal_state = assemblage[goal]

    # Get addition distance
    add_dist = sum(ceil.(max.(goal_state.x - current_state.x, 0.0) / perturbations.epsilon))
    
    # Get deletion distance
    del_dist = sum(ceil.(max.(current_state.x - goal_state.x, 0.0) / perturbations.epsilon))

    # Get temperature distance
    temp_dist = Int64(goal_state.temp != current_state.temp)
    
    return add_dist * costs.add_cost + del_dist * costs.del_cost + temp_dist * costs.temp_cost
end

function net_species_change(current::Int64, 
    goal::Int64, 
    assemblage::Dict{Int64,LVState})
    #= Get net changes between state 1 and 2. 

    Args:
        start: Start state index
        goal: Goal state index
        assemblage: Assemblage dictionary
    
    Return: (tuple)
        net_add
        net_del
        net_temp
    =#
    # Get the LVStates
    current_state = assemblage[current]
    goal_state = assemblage[goal]

    # Get addition distance
    net_add = length(findall(x->x>0,(max.(goal_state.state - current_state.state, 0.0))))
    
    # Get deletion distance
    net_del = length(findall(x->x>0,(max.(current_state.state - goal_state.state, 0.0))))

    # Get temperature distance
    net_temp = Int64(goal_state.temp != current_state.temp)
    
    return (net_add, net_del, net_temp)
end