########################### A* Simulation ###########################
# Running experiment setup -- general
save_data = true  # Saving data to the EcologicalNavigation/data folder
parallelize = true  # Parallelize tasks
procs_num = 4  # Parallel thread number
datetime_string = Dates.format(now(), "e_d_u_Y_HH_MM")
max_candidate_pairs = 10000

# Grid of parameters
add_cost_grid = exp10.(range(-1; stop=1, length=3))
del_cost_grid = exp10.(range(-1; stop=1, length=3))
wait_cost_grid = exp10.(range(-1; stop=1, length=3))
temp_cost_grid = exp10.(range(-1; stop=1, length=3))
epsilon_grid = exp10.(range(-5; stop=-1, length=3))
max_action_grid = [2, typemax(Int64)]

# Experimental data setup:
# Choose among: 
# "Venturelli", "Bucci", "Maynard", "Carrara", 
# "Maynard15-19-23", "Maynard15-17-19-21-23"
experimental_data_set_names = [
	"Venturelli", "Bucci", "Maynard", "Carrara", 
	"Maynard15-19-23", "Maynard15-17-19-21-23"]