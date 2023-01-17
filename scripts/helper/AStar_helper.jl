# Importing all the packages needed for the helper functions
using DataFrames
using Distributions
using LinearAlgebra
using EcologicalNavigation
using Random
using LightGraphs
using SimpleWeightedGraphs
using SparseArrays
using StatsBase
using Printf
using Dates
using IterTools
using SharedArrays
using ProgressMeter
using CSV
using DifferentialEquations

struct ExperimentalData
	#= ExperimentalData struct that contains the setup
    =#
    n_species::Int64
    n_temps::Int64
    data_name::String
	temperature_list::Array{String, 1}
	epsilon::Float64
    data_string::String
end

function load_data(n_species::Int64, 
	n_temps::Int64, 
	data_name::String,
	temperature_list::Vector{String},
	data_directory::String)
	#= Load up experimental or synthetic data.
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
		data_name: Name of dataset
		temperature_list: Temperatures for the dataset (empty string 
				means no specific temperature)
		data_directory: String for the relative/absolute directory
	
	Return:
		A_matrices: A matrices for all temperatures
		r_vectors: r vectors for all temperatures 
		species_names: Names of the species
	=#
	# Create A, r lists
	A_matrices = Matrix{Float64}[]
	r_vectors = Vector{Float64}[]
	species_names = String[]

	# Load experimental data
	# Iterate over all temperatures
	for temp = temperature_list
		# Prepare dataset name
		data_string = data_directory * "/" * data_name
		if temp != ""
			data_string *= "/" * temp
		end

		# Process and push into the list
		A_df = CSV.read(data_string * "/a_matrix.csv", DataFrame)
		r_df = CSV.read(data_string * "/r_vector.csv", DataFrame)
		A = Matrix(A_df)
		r = Vector(r_df[1, :])
		species_names = names(A_df)

		push!(A_matrices, A)
		push!(r_vectors, r)
	end

	return (A_matrices, r_vectors, species_names)
end

function get_data(n_species::Int64, 
	n_temps::Int64, 
	data_string::String,
	A_matrices::Vector{Matrix{Float64}}, 
	r_vectors::Vector{Vector{Float64}};
	show_print_statements = true)
	#= Load up LV parameter settings.
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
		A_matrices: A matrices for all temperatures
		r_vectors: r vectors for all temperatures 

	Keyword Args:
		show_print_statements: If true, shows all print statements
	
	Return:
		lv_params: An LVParamsFull object that contains A, r for all temperatures
	=#
	lv_params = LVParamsFull(A_matrices, r_vectors, 1:n_species)
	
	if show_print_statements
		println("==================================")
		println("Setting up LV system with N = " 
			* string(n_species) * ", T = " * string(n_temps) * "...")
		println("Loading A, r from data: "*data_string)
	end

	if show_print_statements
		println("==================================")
		println("A matrices:")
		for i=1:n_temps
			println("(T = ", i, ")")
			show(stdout, "text/plain", lv_params.A_matrices[i])
			println()
		end
		
		println("\n==================================")
		println("r vectors:")
		for i=1:n_temps
			println("(T = ", i, ")")
			show(stdout, "text/plain", lv_params.r_vectors[i])
			println()
		end
	end

	return lv_params
end

function get_assemblage_transitions(n_species::Int64, 
	n_temps::Int64, 
	perturbations::PerturbationParams,
	costs::TransitionCosts,
	lv_params::LVParamsFull;
	show_print_statements = true,
    ode = ODEParams(Rodas4P(), 1e-6, 1e-3, 1e3))
	#= Generate assemblage and transition
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
        perturbations: Perturbation parameters
        costs: Action costs
		lv_params: An LVParamsFull object that contains 
				A, r for all temperatures
	
	Keyword Args:
		show_print_statements: If true, shows all print statements
        ode: ODE Parameters

	Return:
		assemblage: Assemblage dictionary 
		transitions: Transitions dictionary
		assemblage_df: Assemblage dataframe 
		transitions_df: Transitions dataframe
	=#
	# Generate transition networks
	if show_print_statements
		println("\n==================================")
		println("Generating assemblages and networks...")
	end

	# Generate assemblages
	assemblage = generate_assemblage(
		n_species, n_temps, lv_params;
		show_print_statements = show_print_statements)
	assemblage_df = tabular_assemblage(assemblage)

	# Generate transitions
	transitions_df = generate_lv_transitions(
		n_species, n_temps, assemblage, 
		perturbations, lv_params; 
		remove_loop = true,
		show_print_statements = show_print_statements,
		parallelize = false,
		ode = ode)
	(transitions, _) = dictionary_transitions(transitions_df, costs)

	return (assemblage, transitions, 
		assemblage_df, transitions_df)
end

function get_candidates(assemblage_df::DataFrame;
	show_print_statements = true,
	max_candidate_pairs = typemax(Int64))
	#= Generate candidate states and pairs
		
	Args:
		assemblage_df: Assemblage dataframe
	
	Keyword Args:
		show_print_statements: If true, shows all print statements
		max_candidate_pairs: Maximum candidate pairs

	Return:
		proportion: Proportion of states that are candidates 
				(feasible and stable) among all states
		candidate_states: Vector of all candidate states
		candidate_states_pairs: Vector of all candidate states pairs 
				(i, j) (from i to j)
	=#
	# Test how much improvement there is by testing all the 
	# path finding between the stable population arrangements
	if show_print_statements
		println("==================================")
		println("Finding suitable candidates that are stable "
			* "and feasible and iterating over all path combinations...")
	end

	candidate_states = findall(assemblage_df.candidate)
	proportion = length(candidate_states) * 1.0 / nrow(assemblage_df)
	
	if show_print_statements
		println("Proportion of candidate states: ", 
			proportion)
		println("Candidate states: ", length(candidate_states))
		println("Total states: ", nrow(assemblage_df))
	end

	# Pick all candidate pairs
	candidate_states_pairs = Vector{Int64}[]
	for i = 1:length(candidate_states), j = 1:length(candidate_states)
		if i != j
			push!(
				candidate_states_pairs, 
				[candidate_states[i], candidate_states[j]])
		end
	end

	# Pick a subset of candidate pairs
	if length(candidate_states_pairs) > max_candidate_pairs
		candidate_pairs_indices = StatsBase.sample(
			MersenneTwister(1),
			1:length(candidate_states_pairs), 
			max_candidate_pairs;
			replace = false, ordered = true)
		candidate_states_pairs = candidate_states_pairs[candidate_pairs_indices]
	end

	return (proportion, candidate_states, candidate_states_pairs)
end

function run_a_star(n_species::Int64, 
    n_temps::Int64, 
	candidate_states::Vector{Int64},
    candidate_states_pairs::Vector{Vector{Int64}},
	transitions::Dict{Int64, Vector{LVTransition}},
	assemblage::Dict{Int64,LVState},
    perturbations::PerturbationParams,
	costs::TransitionCosts;
	show_print_statements = true)
	#= Run A* experiment for a given parameter set.
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
        candidate_states: Vector of all candidate states
		candidate_states_pairs: Vector of all candidate states pairs 
				(i, j)
        transitions: Dict{Int64 => Vector{LVTransition}} object that contains 
				all transitions information
		assemblage: Dict{Int64 => LVState} object that contains 
				all assemblage information
        perturbations: Perturbation parameters
        costs: Action costs
		
	Keyword Args:
		show_print_statements: If true, shows all print statements

	Return:
		a_star_df: A DataFrame that contains the
			results of the A* search for going from i to j.
		transition_results: DataFrame of transition statistics
	=#
	# Define the heuristic
	heuristic(x, y) = optimistic_distance(x, y, assemblage, costs)

	# Run the experiments
	if show_print_statements
		results = @showprogress map(states -> solve_a_star(
			states, transitions, heuristic), candidate_states_pairs
			)
	else
		results = map(states -> solve_a_star(
			states, transitions, heuristic), candidate_states_pairs
			)
	end

	# Concatenate results 
	results_table = vcat(results...)

	# Visualize A* results
	a_star_df, transition_results = visualize_a_star_results(n_species, n_temps, 
    	results_table.path, results_table.action, results_table.cost, 
		results_table.idx_pair, results_table.success,
		candidate_states, assemblage, transitions, perturbations, costs)

	return a_star_df, transition_results
end

function solve_a_star(idx_pair::Vector{Int64}, 
	transitions::Dict{Int64, Vector{LVTransition}},
    heuristic::Function)
	#= Function to run a single A* experiment 
		between states i and j.
		
	Args:
		idx_pair: [start, goal] The indices to be tested 
        transitions: Dict{Int64 => Vector{LVTransition}} object that contains 
				all transitions information
        heuristic: Heuristic function that takes in args (index1, index2)
	
	Return: 
		results: Single row DataFrame to be concatenated
	=#
	# Make a data frame for this one thread
	results = DataFrame(
		success = Bool[],
		path = Vector{Int64}[],
		action = Vector{Int64}[], 
		cost = Float64[],
		idx_pair = Vector{Int64}[])

	# Run A*
	start = idx_pair[1]
	goal = idx_pair[2]
	(success, path, action, cost) = a_star_search(start, goal, transitions, heuristic)

	# Get results
	push!(results, [success, path, action, cost, idx_pair])
	
	return results
end

function get_statistics(a_star_df::DataFrame;
	show_print_statements = true)
	#= Calculate statistics for the A* experiment
		
	Args:
		a_star_df: A DataFrame that contains the
			results of the A* search for going from i to j.
        assemblage: Assemblage dictionary
	
	Keyword Args:
		show_print_statements: If true, shows all print statements
	
	Return:
		results_stats: A Dictionary that contains the
			relevant statistics of the results.
	=#
	# Calculate statistics
	a_star_df_filtered = a_star_df[a_star_df.full_path .!= "N/A", :]
	existing_path = nrow(a_star_df_filtered)
	mean_planning = mean(a_star_df_filtered.a_star_cost)
	sd_planning = std(a_star_df_filtered.a_star_cost)
	mean_length = mean(a_star_df_filtered.a_star_ops)
	sd_length = std(a_star_df_filtered.a_star_ops)
	mean_addl_add = mean(a_star_df_filtered.a_star_add - a_star_df_filtered.net_add)
	sd_addl_add = std(a_star_df_filtered.a_star_add - a_star_df_filtered.net_add)
	mean_addl_del = mean(a_star_df_filtered.a_star_del - a_star_df_filtered.net_del)
	sd_addl_del = std(a_star_df_filtered.a_star_del - a_star_df_filtered.net_del)
	mean_addl_temp = mean(a_star_df_filtered.a_star_temp - a_star_df_filtered.net_temp)
	sd_addl_temp = std(a_star_df_filtered.a_star_temp - a_star_df_filtered.net_temp)

	if show_print_statements
		println("==================================")
		println("Printing some interesting statistics...")
		println("- Total # of pairs evaluated: ", nrow(a_star_df))
		println("- Total # of available paths: ", existing_path)
		println("- A* Planning Cost: ", mean_planning, " ± ", sd_planning)
		println("- A* Planning Length: ", mean_length, " ± ", sd_length)
		println("- A* Additional Adds: ", mean_addl_add, " ± ", sd_addl_add)
		println("- A* Additional Dels: ", mean_addl_del, " ± ", sd_addl_del)
		println("- A* Additional Temps: ", mean_addl_temp, " ± ", sd_addl_temp)
	end

	results_stats = Dict(
		"existing_path" => existing_path,
		"mean_planning" => mean_planning,
		"sd_planning" => sd_planning,
		"mean_length" => mean_length,
		"sd_length" => sd_length,
		"mean_addl_add" => mean_addl_add,
		"sd_addl_add" => sd_addl_add,
		"mean_addl_del" => mean_addl_del,
		"sd_addl_del" => sd_addl_del,
		"mean_addl_temp" => mean_addl_temp,
		"sd_addl_temp" => sd_addl_temp,
	)

	return results_stats
end

function save_all_data(n_species::Int64, 
	n_temps::Int64, 
	perturbations::PerturbationParams,
	costs::TransitionCosts,
	lv_params::LVParamsFull,
	data_string::String,
	datetime_string::String,
	species_names::Vector{String},
	temperature_list::Vector{String},
	max_action::Int64,
	assemblage_df::DataFrame,
    transitions_df::DataFrame,
	a_star_results_table::DataFrame,
	a_star_results_stats::Dict,
	transition_results_table::DataFrame,
	candidate_states::Vector{Int64},
	candidate_states_pairs::Vector{Vector{Int64}},
	relative_path::String;
	show_print_statements = true,
	save_assemblage_adjacency = true)
	#= Save the relevant data for the experiment
		
	Args:
		n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
        perturbations: Perturbation parameters
        costs: Action costs
		lv_params: An LVParamsFull object that contains 
				A, r for all temperatures
		data_string: String that indicates the data source
		datetime_string: The string for date-time
		species_names: Vector of strings that show names of species
		temperature_list: Vector of strings of temperature values used
		max_action: Max number of adds and dels
		assemblage_df: Dataframe Assemblage object
        transitions_df: Dataframe of all transitions
		a_star_results_table: A DataFrame that contains the
				results of the A* search for going from i to j.
		a_star_results_stats: A Dictionary that contains the
				relevant statistics of the results.
		transition_results_table: DataFrame of transition statistics
		candidate_states: Vector of all candidate states
		candidate_states_pairs: Vector of all candidate states pairs 
				(i, j) (from i to j)
		relative_path: Relative path to store the data into
	
	Keyword Args:
		show_print_statements: If true, shows all print statements
		save_assemblage_adjacency: If true, save the table for assemblage
				and transtions. Flag it so that we only do this once
				due to multiple grid sweeps over different costs,
				which doesn't change the assemblage and transition edges themselves.

	=#
	# Formatting output file name
	configure_string = "eps" * replace(string(perturbations.epsilon), "." => "-")
	if max_action == typemax(Int64)
		configure_string *= "_uncapped"
	else
		configure_string *= "_capped_" * string(max_action)
	end

	coststring = configure_string * (
		"_add" * replace(string(costs.add_cost), "." => "-") *
		"_del" * replace(string(costs.del_cost), "." => "-") *
		"_temp" * replace(string(costs.temp_cost), "." => "-") *
		"_wait" * replace(string(costs.wait_cost), "." => "-"))

	# Setting up file directories
	fname_directory_string = joinpath(
		relative_path, 
		data_string,
		datetime_string)
	isdir(fname_directory_string) || mkpath(fname_directory_string)
	
	# Save A star and transitions
	CSV.write(fname_directory_string * "/" * coststring * "_a_star.csv", a_star_results_table)
	CSV.write(fname_directory_string * "/" * coststring * "_transitions.csv", transition_results_table)

	if save_assemblage_adjacency
		transitions_df_print = select(transitions_df, Not(:cost))
		CSV.write(fname_directory_string * "/" * configure_string * "_assemblage.csv", assemblage_df)
		CSV.write(fname_directory_string * "/" * configure_string * "_adjacency_df.csv", transitions_df_print)
	end
	
	if show_print_statements
		println(stderr,"\n* Saving experimental results to ", fname_directory_string)
	end

	# Writing output files
	open(fname_directory_string * "/" * coststring * "_summary.txt", "a") do io
		println(io, "==================================")
		println(io, "Setup:")
		println(io, "N: ", n_species)
		println(io, "T: ", n_temps)
		println(io, "add: ", costs.add_cost)
		println(io, "del: ", costs.del_cost)
		println(io, "wait: ", costs.wait_cost)
		println(io, "temp: ", costs.temp_cost)
		println(io, "integration_time_multiplier: ", perturbations.time)
		println(io, "epsilon: ", perturbations.epsilon)
		println(io, "tolerance: ", perturbations.tolerance)
		println(io, "Data set: ", data_string)
		println(io, "Species names: ", species_names)
		println(io, 
			"Species label: ", 
			[x for x = 1:length(species_names)])
		println(io, "Temperature names: ", temperature_list)
		println(io, 
			"Temperature label: ", 
			[x for x = 1:length(temperature_list)])
		
		println(io)

		println(io, "==================================")
		println(io, "Statistics:")
		println(io, 
			"Proportion of candidate states: ", 
			length(candidate_states)*1.0/nrow(assemblage_df))
		println(io, "Candidate states: ", length(candidate_states))
		println(io, "Total states: ", nrow(assemblage_df))
		println(io, "Total natural transitions (may include duplicates): ", nrow(transitions_df))
		println(io, "Total candidate pairs evaluated: ", 
			length(candidate_states_pairs))
		println(io, "Total viable paths: ", a_star_results_stats["existing_path"])
		replace!(a_star_results_table.a_star_ops, missing => -1)
		println(io, "Total viable paths (non-trivial length): ", 
			nrow(a_star_results_table[a_star_results_table.a_star_ops.>1,:]))
		println(io, "Total viable path proportions (viable / # candidate pairs): ", 
			a_star_results_stats["existing_path"]*1.0/length(candidate_states_pairs))
		println(io, "Total non-trivial paths proportions (non-trivial / viable): ", 
			nrow(a_star_results_table[a_star_results_table.a_star_ops.>1,:])*1.0/a_star_results_stats["existing_path"]) 
		println(io)
		println(io, "Planning Cost - Mean: ", a_star_results_stats["mean_planning"])
		println(io, "Planning Cost - SD: ", a_star_results_stats["sd_planning"])
		println(io, "Planning Length - Mean: ", a_star_results_stats["mean_length"])
		println(io, "Planning Length - SD: ", a_star_results_stats["sd_length"])
		println(io)
		println(io, "Additional Adds - Mean: ", a_star_results_stats["mean_addl_add"])
		println(io, "Additional Adds - SD: ", a_star_results_stats["sd_addl_add"])
		println(io, "Additional Dels - Mean: ", a_star_results_stats["mean_addl_del"])
		println(io, "Additional Dels - SD: ", a_star_results_stats["sd_addl_del"])
		println(io, "Additional Temps - Mean: ", a_star_results_stats["mean_addl_temp"])
		println(io, "Additional Temps - SD: ", a_star_results_stats["sd_addl_temp"])
		println(io)

		println(io, "==================================")
		println(io, "A matrix - Entry statistics:")
		for i=1:n_temps
			println(io,"(T = ", i, ")")
			println(io, "A matrix - Mean: ", mean(lv_params.A_matrices[i]))
			println(io, "A matrix - SD: ", std(lv_params.A_matrices[i]))
		end
		println(io, "r vector - Entry statistics:")
		for i=1:n_temps
			println(io,"(T = ", i, ")")
			println(io, "r vector - Mean: ", mean(lv_params.r_vectors[i]))
			println(io, "r vector - SD: ", std(lv_params.r_vectors[i]))
		end

		println(io, "==================================")
		println(io, "A matrices:")
		for i=1:n_temps
			println(io,"(T = ", i, ")")
			show(io, "text/plain", lv_params.A_matrices[i])
			println(io)
		end

		println(io, "\n==================================")
		println(io, "r vectors:")
		for i=1:n_temps
			println(io, "(T = ", i, ")")
			show(io, "text/plain", lv_params.r_vectors[i])
			println(io)
		end
	end

	if show_print_statements
		println(stderr, "* Saving summary to ", fname_directory_string)
	end
end

function run_full_experiments(n_species::Int64,
    n_temps::Int64,
    data_name::String,
	temperature_list::Vector{String},
	epsilon::Float64,
    data_string::String,
	datetime_string::String,
	add_cost_grid::Vector{Float64},
	del_cost_grid::Vector{Float64},
	wait_cost_grid::Vector{Float64},
	temp_cost_grid::Vector{Float64},
	max_action_grid::Vector{Int64};
	save_data = false,
	relative_path = "",
	max_candidate_pairs = typemax(Int64))
	#= Run the full A* experiment.
		
	Args:
		n_species: Number of species in the system
		n_temps: Number of temperature settings in the system
		data_name: Name of dataset
		temperature_list: Temperatures for the dataset (empty string 
				means no specific temperature)
		epsilon: Integration parameter
		data_string: The raw string for dataset
		datetime_string: The string for date-time
		add_cost_grid: Vector of grid values for add
		del_cost_grid: Vector of grid values for del
		wait_cost_grid: Vector of grid values for wait
		temp_cost_grid: Vector of grid values for temp
		max_action_grid: Vector of grid values for max_action
	
	Keyword Args:
		save_data: Boolean for whether to save data files
		relative_path: Relative path to save files to
		max_candidate_pairs: Maximum candidate pairs

	Return:
		None -- saved to all files
	=#
	# Set parameters for the dataset
    perturbations = PerturbationParams(5, epsilon, 5e-2)
	ode = ODEParams(Rodas4P(), 1e-6, 1e-6, 1e3)

    # Load up system parameters
    (A_matrices, r_vectors, species_names) = load_data(
            n_species, n_temps, data_name, 
            temperature_list, dataset_string)
    
    # Setup the environment
    lv_params = get_data(
        n_species, n_temps, data_string, A_matrices, r_vectors;
        show_print_statements = true)
    temp_costs = TransitionCosts(1.0, 1.0, 1.0, 1.0)
    (assemblage, _, assemblage_df, transitions_df) = get_assemblage_transitions(
        n_species, n_temps, perturbations, temp_costs, lv_params; 
        show_print_statements = true,
        ode = ode)
    (proportion, candidate_states, candidate_states_pairs) = get_candidates(
        assemblage_df;
        show_print_statements = true,
        max_candidate_pairs = max_candidate_pairs)
    
    # Get all costs
    for max_action in max_action_grid
		first_loop = true
	    for addc in add_cost_grid, 
	    	delc in del_cost_grid, 
		    waitc in wait_cost_grid, 
		    tempc in temp_cost_grid

	        # Regenerate transition dictionary with new costs
	        costs = TransitionCosts(addc, delc, waitc, tempc)
	        (transitions, transitions_df_filtered) = dictionary_transitions(
	        	transitions_df, costs; max_action = max_action)

	        # Run the experiments
	        a_star_results_table, transition_results_table = run_a_star(
	            n_species, n_temps, candidate_states,
	            candidate_states_pairs, transitions, assemblage, 
	            perturbations, costs;
	            show_print_statements = true)
	        a_star_results_stats = get_statistics(
	            a_star_results_table;
	            show_print_statements = true)
	        
	        # Save Data
	        if save_data
	            # Make directory if it does not exist
	            isdir(relative_path) || mkpath(relative_path)

	            # Save all the data
	            save_all_data(n_species, n_temps, perturbations, costs, lv_params, 
					data_string, datetime_string, species_names, temperature_list, max_action,
	                assemblage_df, transitions_df_filtered, a_star_results_table, a_star_results_stats, 
	                transition_results_table, candidate_states, candidate_states_pairs, relative_path; 
	                show_print_statements = true,
					save_assemblage_adjacency = first_loop)
				first_loop = false
	        end
		end
	end
end

run_full_experiments(experiment_params::ExperimentalData, 
	relative_path::String,
	datetime_string::String,
	add_cost_grid::Vector{Float64},
	del_cost_grid::Vector{Float64},
	wait_cost_grid::Vector{Float64},
	temp_cost_grid::Vector{Float64},
	max_action_grid::Vector{Int64},
	save_data::Bool,
	max_candidate_pairs::Int64) = run_full_experiments(
		experiment_params.n_species,
		experiment_params.n_temps,
		experiment_params.data_name,
		experiment_params.temperature_list,
		experiment_params.epsilon,
		experiment_params.data_string,
		datetime_string,
		add_cost_grid,
		del_cost_grid,
		wait_cost_grid,
		temp_cost_grid,
		max_action_grid;
		save_data = save_data,
		relative_path = relative_path,
		max_candidate_pairs = max_candidate_pairs)