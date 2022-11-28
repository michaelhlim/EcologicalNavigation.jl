########################### Setup for the experiments ###########################
# Loading packages
using ClusterManagers
using Distributed

# Set up directory -- SLURM uses different directory access
if "SLURM_JOBID" in keys(ENV)
	directory_string = "./helper"
	dataset_string = "./data/dataset"
    relative_path_string = "./data/results"
else 
	directory_string =  "./helper" 
	dataset_string = "./data/dataset"
    relative_path_string = "./data/results"
end

# Load configs & helpers
include(directory_string * "/AStar_helper.jl")
include(directory_string * "/experiment_config.jl")

# Loading packages for every compute node to parallelize
if parallelize
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()

	worker_ids = if "SLURM_JOBID" in keys(ENV)
		wids = ClusterManagers.addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
		wids
	else
		addprocs(procs_num)
	end  

	Distributed.@everywhere begin
		using Pkg
		Pkg.activate(".")
		Pkg.instantiate()
	end

	Distributed.@everywhere begin
		# Set up directory -- Need to do this again...
		if "SLURM_JOBID" in keys(ENV)
			Distributed.@everywhere directory_string = "./helper"
			Distributed.@everywhere dataset_string = "./data/dataset"
            Distributed.@everywhere relative_path_string = "./data/results"
		else 
			Distributed.@everywhere directory_string =  "./helper" 
			Distributed.@everywhere dataset_string = "./data/dataset"
            Distributed.@everywhere relative_path_string = "./data/results"
		end

		# Load helper functions
		include(directory_string * "/AStar_helper.jl")
	end
	@info "Started $(Distributed.nworkers()) workers..."
	@info "Precompiling simulation code..."

end


########################### Simulation ###########################
# Run A* for all dataset
try
	include(directory_string * "/load_data.jl") 
    if parallelize
        pmap(col -> run_full_experiments(col, relative_path_string, datetime_string, 
			add_cost_grid, del_cost_grid, wait_cost_grid, temp_cost_grid, max_action_grid, 
			save_data, max_candidate_pairs), 
			experiment_params)
    else
        map(col -> run_full_experiments(col, relative_path_string, datetime_string, 
			add_cost_grid, del_cost_grid, wait_cost_grid, temp_cost_grid, max_action_grid, 
			save_data, max_candidate_pairs), 
			experiment_params)
    end
finally
	if parallelize
		Distributed.rmprocs(worker_ids)
	end
end