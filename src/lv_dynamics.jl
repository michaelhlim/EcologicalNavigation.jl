#= Lotka-Volterra script that contains all the generation and configuration for assemblages.
=#
struct LVParams
    #= LVParams struct that contains the A, r matrices 
        for a given system (fixed T).
        
    Args:
        A_matrix: A matrix
        r_vector: r vector
        collabels: Label of each columns, usually useful for determining 
                which species are actually present when dealing with 
                subsystem of A, r
    =#
    A_matrix::Array{Float64, 2}
    r_vector::Array{Float64, 1}
    collabels::Array{Float64, 1}
end

struct LVParamsFull
    #= LVParamsFull struct that contains the A, r matrices 
        for a given system (multiple T).
        
    Args:
        A_matrices: Array of A matrices for each temperature
        r_vectors: Array of r vectors for each temperature
        collabels: Label of each columns, usually useful for determining 
                which species are actually present when dealing with 
                subsystem of A, r
    =#
    A_matrices::Array{Array{Float64, 2}, 1}
    r_vectors::Array{Array{Float64, 1}, 1}
    collabels::Array{Float64, 1}
end

struct LVState
    #= LVState struct that contains all the information about a state configuration.
        
    Args:
        state: A binary vector of whether a species is existent
        x: A real vector of the actual abundance of each species
        temp: Temperature variable discretized into T bins
        str: String representation of the state
        stable: Indicator for whether the state is stable
        feasible: Indicator for whether the state is feasible
        candidate: Indicator for whether the state is both stable and feasible
        tau: The time constant, if exists
        richness: Total number of unique existing species
        abundance_mean: Mean species abundance for existing species
        abundance_sd: Stdev of species abundance for existing species
    =#
    state::Array{Int64, 1}
    x::Array{Float64, 1}
    temp::Int64
    str::String
    stable::Bool
    feasible::Bool
    candidate::Bool
    tau::Float64
    richness::Int64
    abundance_mean::Float64
    abundance_sd::Float64
end

struct LVAction
    #= LVAction struct that contains the state and temperature changes.
        
    Args:
        ds: State change: -1 is deletion, 0 is no action, and 1 is addition
        dt: Temperature change to a desired temperature
        str: String representation of the action
    =#
    ds::Array{Int64, 1}
    dt::Int64
    str::String
end

struct LVTransition
    #= LVTransition struct that contains a transition between different states.
        
    Args:
        s: Starting state index
        a: Action index
        f: Finishing state index
        cost: Cost of the action
        a_str: Action string
    =#
    s::Int64
    a::Int64
    f::Int64
    cost::Float64
    a_str::String
end

struct FeasibilityStruct
    #= FeasibilityStruct struct that contains the information used to 
        calculate feasible of a state.
        
    Args:
        feasible: Boolean for indicating whether the struct is feasible
        new_state: Indicator vector for the species have nonzero survival
        equilibrium_state: The actual equilibrium vector of the survival state
        abundance_mean: Mean species abundance for existing species
        abundance_sd: Stdev of species abundance for existing species
    =#
    feasible::Bool
    new_state::Array{Int64, 1}
    equilibrium_state::Array{Float64, 1}
    abundance_mean::Float64
    abundance_sd::Float64
end

struct StabilityStruct
    #= StabilityStruct struct that contains the information used to 
        calculate stable of a state.
        
    Args:
        stable: Boolean for indicating whether the struct is stable
        lambda: Eigenvalue of D(x^*)*A, where x^* is the equilibrium point.
                x^* = -A^-1*r = 1/T \int_0^T x(t) dt
        tau: The time constant, which is inverse of the biggest eigenvalue 
                if exists
    =#
    stable::Bool
    lambda::Array{Float64, 1}
    tau::Float64
end

struct PerturbationParams
    #= PerturbationParams struct that contains parameters for 
        transitions with perturbations
        
    Args:
        time: Factor to multiply with max_tau to get max time for ODE solver range (min time is 0)
        epsilon: Amount of addition and deletion
        tolerance: Tolerance for comparing integration vs. analytical equilibrium
    =#
    time::Float64
    epsilon::Float64
    tolerance::Float64
end

struct TransitionCosts
    #= TransitionCosts struct that contains the costs for a given system.
        
    Args:
        add_cost: Cost of addition
        del_cost: Cost of deletion
        wait_cost: Cost of wait
        temp_cost: Cost of temperature change
    =#
    add_cost::Float64
    del_cost::Float64
    wait_cost::Float64
    temp_cost::Float64
end

struct ODEParams
    #= ODEParams struct that contains the settings for ODE.
        
    Args:
        alg: ODE Solver algorithm
        abstol: Absolute tolerance
        reltol: Relative tolerance
        maxiters: Max iterations
    =#
    alg::Any
    abstol::Float64
    reltol::Float64
    maxiters::Float64
end

function generate_states(n_species::Int64, n_temps::Int64)
    #= Create all possible 2^n_species * n_temps state vectors.
        
    Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
    
    Return: (tuple)
        state_vecs: An Array of all possible states
        state_temps: An Array of all possible temps
    =#
    # Get all states with binary operation
    state_vecs = Array{Int64,1}[]
    for i = 0:2^(n_species)-1
        s = bitstring(i)
        s = reverse(s[end-n_species+1:end])
        s_vec = [parse(Int64, ss) for ss in split(s, "")]
        push!(state_vecs, s_vec)
    end
    state_vecs = repeat(state_vecs, n_temps)
    
    # Get all the temperatures
    state_temps = repeat(1:n_temps, inner=[2^n_species])

    return (state_vecs, state_temps)
end

function determine_feasible(lv_params::LVParams)
    #= Determine feasible of the (sub)parameter set.
    Feasibility is calculated by seeing if the resulting
    equilibrium state x^* = -A^-1*r has all positive values.
    Also returns the equilibrium state itself, as well as
    the abundance mean and sd.
        
    Args:
        lv_params: LVParams object of (sub)parameters

    Return:
        FeasibilityStruct that contains the feasible 
        information of the parameter set
    =#
    # Catch if the r vector is null
    if length(lv_params.r_vector) > 0
        # New abundance vector via LV calculation
        x_vector = -1.0 * inv(lv_params.A_matrix) * lv_params.r_vector

        # Check if the solution has all positive values
        feasible = !any(x->x<=0, x_vector)
        
        # Get the new state where the species 
        # have positive abundance
        new_state_idx = findall(x -> x>0, x_vector)
		new_state = lv_params.collabels[new_state_idx]
        equilibrium_state = x_vector

        # Get mean and sd of residents
        abundance_mean = mean(x_vector[new_state_idx])
        abundance_sd = std(x_vector[new_state_idx])
    else
        feasible = true
        new_state = Int64[]
        equilibrium_state = Float64[]
        abundance_mean = NaN
        abundance_sd = NaN
    end

    return FeasibilityStruct(
        feasible, 
        new_state, 
        equilibrium_state,
        abundance_mean, 
        abundance_sd)
end

function determine_stable(lv_params::LVParams)
    #= Determine stable of the (sub)parameter set.
    Stability is calculated by seeing if the eigenvalues of
    D(x^*)*A are all negative, with equilibrium state
    x^* = -A^-1*r and D() the diagonalization operator.
    Also returns the tau, time constant.
        
    Args:
        lv_params: LVParams object of (sub)parameters

    Return:
        Stability struct that contains the stable 
        information of the parameter set
    =#
    if length(lv_params.r_vector) > 0
        # Equilibrium state
        x_vector = -1.0 * inv(lv_params.A_matrix) * lv_params.r_vector

        # Calculate eigenvalues and see if 
        # max(Re(lambda)) < 0 for all lambda's
        lambda = eigvals(Diagonal(x_vector) * lv_params.A_matrix)
        if typeof(lambda[1]) == Complex{Float64}
            lambda = [z.re for z in lambda]
        end
        stable = !any(x->x>=0, lambda)

        # Time constant is inverse of 
        # the biggest eigenvalue, if exists
        if stable
            tau = -1/maximum(lambda) 
        else
            tau = 0
        end
    else
        lambda = Float64[]
        stable = true
        tau = NaN
    end

    return StabilityStruct(stable, lambda, tau)
end

function generate_assemblage(n_species::Int64, 
	n_temps::Int64, 
	lv_params::LVParamsFull;
	show_print_statements = true)
	#= Generate full assemblage for a given parameter configuration.
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
		lv_params: An LVParamsFull object that contains 
				A, r for all temperatures
	
	Keyword Args:
		show_print_statements: If true, shows all print statements
	
	Return:
		assemblage: Dict{Int64 => LVState} object that contains 
				all assemblage information
	=#
    # Get all the states
    (state_vecs, state_temps) =  generate_states(n_species, n_temps)
    n_assemblage = (2^n_species) * n_temps
    assemblage = Dict{Int64, LVState}()

	if show_print_statements
		println(("Generating assemblage with N = " 
				* string(n_species) * " and T = " * string(n_temps) * "...")
		)
	end
	prog_bar = Progress(n_assemblage; 
		showspeed = true,
		enabled = show_print_statements
		)
	
	for i = 1:n_assemblage # Skip no species assemblage
        state = state_vecs[i]
        temp = state_temps[i]
        str_representation = state_idx_to_str(n_species, n_temps, i)
		lv_params_this = _generate_subparameters(state, temp, lv_params)
        
    	# Determine stable and feasible
	    f = determine_feasible(lv_params_this)
	    s = determine_stable(lv_params_this)
        equilibrium_state = zeros(Float64, n_species)
        if f.feasible
	    	equilibrium_state[f.new_state] = f.equilibrium_state
	    end
	    
        # Insert into the dictionary
        assemblage[i] = LVState(
            state,  # state
            equilibrium_state,  # x
            temp,  # temp
            str_representation,  # str
            s.stable,  # stable
            f.feasible,  # feasible
            s.stable & f.feasible,  # candidate
            s.tau,  # tau
            sum(state),  # richness
            f.abundance_mean,  # abundance_mean
            f.abundance_sd,  # abundance_sd
        )

		next!(prog_bar)
	end
    
	return assemblage
end

function tabular_assemblage(assemblage::Dict{Int64,LVState})
    #= Return tabular form of assemblage for better interpretation
		
	Args:
		assemblage: Dict{Int64 => LVState} object that contains 
				all assemblage information

	Return:
		assemblage_df: Dataframe object that contains 
				all assemblage information
	=#
    assemblage_df = DataFrame(
        state = String[],
        x = String[],
        temp = Int64[],
        str = String[],
        stable = Bool[],
        feasible = Bool[],
        candidate = Bool[],
        tau = Float64[],
        richness = Int64[],
        abundance_mean = Float64[],
        abundance_sd = Float64[],
    )
    
    # Iterate over the dictionary in order of index
    for s in sort(collect(assemblage), by=x->x[1])
        push!(assemblage_df,[
            replace(string(s.second.state), "," => " "),
            replace(string(s.second.x), "," => " "),
            s.second.temp,
            s.second.str,
            s.second.stable,
            s.second.feasible,
            s.second.candidate,
            s.second.tau,
            s.second.richness,
            s.second.abundance_mean,
            s.second.abundance_sd,
        ])
    end

    return assemblage_df
end

function valid_delta_action(state::LVState,
    action::LVAction,
    perturbations::PerturbationParams; 
	remove_loop = true)
    #= Augment state with delta action if possible, 
        otherwise return false
        
    Args:
        state: LVState variable
        action: LVAction variable
        perturbations: Perturbation parameters
    
    Keyword Args:
		remove_loop: Remove any self-moves
    
    Return: (tuple)
        new_s: State after augmenting with action
        new_t: Temperature after augmenting with action
        num_add: Number of additions
        num_del: Number of deletions
        num_temp: Number of temperature changes 
        action_string: String of performed actions
        valid: Boolean for if action is valid
    =#
    # Check if deletions are valid
    if ~issubset(findall(x -> x<0, action.ds), findall(x -> x>0, state.x))
        return (state.x, state.temp, 0, 0, 0, "", false)
    end

    # Check if ds is zero and there is no temp change
    if !any(x -> x!=0, action.ds) && (state.temp == action.dt)
        return (state.x, state.temp, 0, 0, 0, "", false)
    end

    # Set up tracking variables
    num_add = 0
    num_del = 0
    num_temp = 0

    # Apply the actions
    s_vec = copy(state.x)
    
    # Additions
    adds = findall(x -> x>0, action.ds)
    s_vec[adds] .+= perturbations.epsilon
    num_add += length(adds)
    
    # Deletions
    dels = findall(x -> x<0, action.ds)
    s_vec[dels] .-= perturbations.epsilon
    s_vec = max.(s_vec, 0.0)
    num_del += length(dels)

    # Check if temperature changed
    temp_string = ""
    if state.temp != action.dt
        temp_string = "(*" * string(action.dt) * ")"
        num_temp += 1
    end

    # Process all strings
    action_string = action.str * temp_string

    return (s_vec, action.dt, num_add, num_del, num_temp, action_string, true)
end

function single_lv_transition(delta_state::Vector{Float64}, 
    time::Float64,
	A_matrix::Matrix{Float64},
    r_vector::Vector{Float64};
    ode = ODEParams(Rodas4P(), 1e-6, 1e-3, 1e3))
    #= Function that generates a single transition from the LVState.

	Args:
		delta_state: A transient state vector
        time: Max time for ODE solver range (min time is 0)
    
    Keyword Args:
        ode: ODE Parameters

	Returns: (tuple)
		new_state: New state vector after the transition
        success: Boolean for success or failure to solve
	=#

    # Get the differential equation
    u0 = delta_state
    time_span = (0.0, time)
    lv_diffeq(u, p, t) = Diagonal(u) * (r_vector + A_matrix * u)

    # println(delta_state, A_matrix, r_vector, time)
    # Solve the differential equation
    problem = ODEProblem(lv_diffeq, u0, time_span)
    solution = solve(problem, alg=ode.alg, 
        abstol=ode.abstol, reltol=ode.reltol, 
        maxiters=ode.maxiters; 
        verbose=false)

    # Return Variables
    new_state = solution.u[end]
    success = (solution.retcode == :Success)

    return (new_state, success)
end

function generate_lv_transitions(n_species::Int64, 
	n_temps::Int64, 
	assemblage::Dict{Int64,LVState}, 
    perturbations::PerturbationParams,
	lv_params::LVParamsFull; 
	remove_loop = true,
	show_print_statements = true,
    parallelize = false,
    ode = ODEParams(Rodas4P(), 1e-6, 1e-3, 1e3))
	#= Function that generates the transition network given an assemblage.
        Parallelized version.

	Args:
		n_species: Number of species (N)
		n_temps: Number of temperature settings (T)
		assemblage: Assemblage dictionary
        perturbations: Perturbation parameters
		lv_params: Parameters, which contains A, r matrices
	
	Keyword Args:
		remove_loop: Remove any self-moves
		show_print_statements: If true, shows all print statements
        parallelize: Decide to run pmap vs. map
        ode: ODE Parameters

	Returns:
		transitions_df: DataFrame, representing all possible AStar edges
	=#
	# Get appropriate constants
    n_actions = ((3^n_species) * n_temps)
    max_tau = 0.0
    for (_, v) in assemblage
        if v.candidate && !isnan(v.tau)
            max_tau = max(max_tau, v.tau)
        end
    end

    if show_print_statements
        println(
            "Generating single species connections for A* graph with N = "
            * string(n_species) * ", T = " * string(n_temps) 
            * ".")
    end

    # Getting the assemblage to create KD tree for each temperature
    # (We shouldn't need to create it for each temperature, but just to be safe)
    assemblage_filtered = Dict{Int64,Vector{Vector{Float64}}}()
    for i = 1:n_temps
        assemblage_filtered[i] = Vector{Float64}[]
    end

    kd_idx = Dict{Tuple{Int64,Int64},Int64}()
    count = ones(Int64, n_temps)
    for i = 1:length(assemblage)
        state = assemblage[i]
        if state.candidate
            push!(assemblage_filtered[state.temp], state.x)
            kd_idx[(state.temp, count[state.temp])] = i
            count[state.temp] += 1
        end
    end
    kd_tree = Dict{Int64, KDTree}()
    for i = 1:n_temps
        kd_tree[i] = KDTree(reduce(hcat, assemblage_filtered[i]))
    end

    # Make a data frame for this one parallel thread
	transitions_df = DataFrame(
		s = Int64[],
        a = Int64[],
        f = Int64[],
        add = Int64[], 
        del = Int64[], 
        temp = Int64[],
        a_str = String[])

    @showprogress for i = 1:length(assemblage), j = 1:n_actions
        state = assemblage[i]
        action = action_idx_to_lv_action(n_species, n_temps, j)
        # Find all possible transitions for all existing entries in the assemblage
        if state.candidate
            # Check if the actions are valid
            (delta_state, delta_temp, add, del, temp, action_string, valid) = valid_delta_action(
                state, action, perturbations; 
                remove_loop = remove_loop)
            if valid
                existent = findall(x->x>0, delta_state)
                if length(existent) > 0
                    lv_sub = _generate_subparameters(delta_state, delta_temp, lv_params)

                    # Obtain new state from LV dynamics integration
                    (new_existent_state, success) = single_lv_transition(
                        delta_state[existent], perturbations.time * max_tau, 
                        lv_sub.A_matrix, lv_sub.r_vector; ode = ode)
                    new_state = zeros(length(delta_state))
                    new_state[existent] = new_existent_state
                    
                    if success
                        # Get the closest equilibrium state candidate
                        idx, dist = nn(kd_tree[action.dt], new_state)
                        to_idx = kd_idx[(action.dt, idx)]
                        candidate_state = assemblage[to_idx]

                        # Check that the equilibrium state matches the new state
                        if dist <= LinearAlgebra.norm(candidate_state.x) * perturbations.tolerance
                            push!(transitions_df, [i, j, to_idx, add, del, temp, action_string])
                        end
                    end
                end
            end
        end
    end

    # Get start and final state strings for better interpretability
    transitions_df.start = map(
        col -> state_idx_to_str(n_species, n_temps, col), 
        transitions_df[!, :s])
    transitions_df.final = map(
        col -> state_idx_to_str(n_species, n_temps, col), 
        transitions_df[!, :f])

    # Remove any self-loops as a result of transition
    if remove_loop
        transitions_df = transitions_df[transitions_df.s .!= transitions_df.f, :]
    end

	return transitions_df
end

function _apply_cost(transitions_df::DataFrame,
    costs::TransitionCosts)
    #= Apply the costs of transitioning for a specific set of parameters
        and returns a copy of the original transitions dataframe
		
	Args:
		transitions_df: Dataframe object that contains 
				all transitions information
        costs: Action costs

	Return:
		transitions_df_cost: Dataframe object that contains 
				all transitions information with cost 
	=#
    transitions_df_cost = copy(transitions_df)
    transitions_df_cost.cost = (
        transitions_df_cost.add * costs.add_cost +
        transitions_df_cost.del * costs.del_cost +
        transitions_df_cost.temp * costs.temp_cost) .+
        costs.wait_cost
    
    return transitions_df_cost
end

function dictionary_transitions(transitions_df::DataFrame,
    costs::TransitionCosts;
    max_action = typemax(Int64))
    #= Return dictionary form of transitions for better computation.
        This will return an adjacency list from the start node
		
	Args:
		transitions_df: Dataframe object that contains 
				all transitions information
        costs: Action costs

    Keyword Args:
        max_action: Maximum amount of adds and dels

	Return:
		transitions: Dict{Int64 => Vector{LVTransition}} object that contains 
				all transitions information
        transitions_df_cost: Filtered transitions dataframe with costs
	=#
    # Apply the costs of transitioning for a specific set of parameters
    transitions_df_cost = _apply_cost(transitions_df, costs)

    # Cap the max number of actions if desired
    transitions_df_cost = filter(
        row -> (row.add <= max_action && row.del <= max_action), 
        transitions_df_cost)
    
    # Prune to get the lowest cost action for a start -> finish
    transitions_df_cost = combine(
        df -> df[argmin(df.cost), :], 
        groupby(transitions_df_cost, [:s, :f]))
    sort!(transitions_df_cost, [:s, :f])

    # Transition network, which is an adjacency list of 
    # Dict{LVState index => List of LVTransition's}
	transitions = Dict{Int64, Vector{LVTransition}}()

    # Make it a dictionary
    start_nodes = unique(transitions_df_cost.s)
    for node in start_nodes
        results_from_node = transitions_df_cost[transitions_df_cost.s .== node, :]
        for row in eachrow(results_from_node)
            transition = LVTransition(
                row.s,
                row.a,
                row.f,
                row.cost,
                row.a_str)
            if haskey(transitions, node)
                push!(transitions[node], transition)
            else
                transitions[node] = [transition]
            end
        end
    end

    return (transitions, transitions_df_cost)
end