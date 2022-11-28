directory_string =  "./helper" 
dataset_string = "./data/dataset"
relative_path_string = "./data/results"

# Load configs & helpers
include(directory_string * "/AStar_helper.jl")
include(directory_string * "/experiment_config.jl")
include(directory_string * "/load_data.jl") 
run_full_experiments(experiment_params[1], relative_path_string, datetime_string, 
	add_cost_grid, del_cost_grid, wait_cost_grid, temp_cost_grid, max_action_grid, 
	save_data, max_candidate_pairs)
