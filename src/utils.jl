#= Utility functions and structs
=#
function _generate_subparameters(state::Vector{Int64}, 
    temperature::Int64, 
    params::LVParamsFull)
	#= Calculate the A, r subparameters according to 
	which states are present in the 'state' vector.
		
	Args:
		state: Vector of species occupancy
        temperature: Temperature index
		params: LVParamsFull parameters set containing LV parameters
				for all temperatures

	Return:
		params: LVParams object that contains subparameters 
				(submatrix, subvector) only for species 
				that exist in the state
	=#
	# Find which species are present in that row
	species_idx = findall(x -> x>0, state)

	# Pick subset of parameters 
	# (assuming that the parameters don't change when subsetting)
	A_subset = [params.A_matrices[temperature][i,j] for i = species_idx, j = species_idx]
	r_subset = params.r_vectors[temperature][species_idx]
  	
  	# Return params
	params = LVParams(
		A_subset, 
		r_subset, 
		species_idx)

	return params
end

function _generate_subparameters(state::Vector{Float64}, 
    temperature::Int64, 
    params::LVParamsFull)
	#= Override
	=#
	# Find which species are present in that row
	species_idx = findall(x -> x>0, state)

	# Pick subset of parameters 
	# (assuming that the parameters don't change when subsetting)
	A_subset = [params.A_matrices[temperature][i,j] for i = species_idx, j = species_idx]
	r_subset = params.r_vectors[temperature][species_idx]
  	
  	# Return params
	params = LVParams(
		A_subset, 
		r_subset, 
		species_idx)

	return params
end

function split_idx(n_species::Int64, n_temps::Int64, idx::Int64)
    #= Split the raw index into a tuple of (species, temp).
        
    Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
        idx: Raw index that varies from 1 to (2^n_species) * n_temps
    
    Return:
        A tuple (idx_species, idx_temp)
        idx_species: Index for species space
        idx_temp: Index for temperature space
    =#
    idx_species = mod(idx-1, 2^n_species) + 1
    idx_temp = Int64(div(idx-1, 2^n_species)) + 1

    return (idx_species, idx_temp)
end

function state_str_to_idx(n_species::Int64, 
    str::AbstractString)
    #= Get state index from string (single temperature).
        
    Args:
        n_species: Number of species in the system
        str: Raw string of state
    
    Return:
        idx: Raw index of state
    =#
    # Empty state is index 1
    if length(str) == 0
        return 1
    end

    # Parse the string
    state_str = parse.(Int64, split(str, "*"))
    idx = 1
    for i in state_str
        idx += 2^(i-1)
    end

    return idx
end

function state_str_to_idx(n_species::Int64, 
    n_temps::Int64, 
    str::AbstractString)
    #= Get state index from string (multi temperature).
        
    Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
        str: Raw string of state
    
    Return:
        idx: Raw index of state
    =#
    str_split = split(str, "|")
    str_species = str_split[1]
    idx_temp = parse(Int64, str_split[2])

    # Get index from sub function
    idx_species = state_str_to_idx(n_species, str_species)

    # Augment by temperature
    idx = idx_species + (idx_temp-1)*2^n_species

    return idx
end

function state_idx_to_vec(n_species::Int64, idx::Int64)
    #= Get state vector from index (single temperature).
        
    Args:
        n_species: Number of species in the system
        idx: Raw index of state
    
    Return:
        s_vec: Vector of state
    =#
    # Remember that indexing starts at 1
    s = bitstring(idx-1)
    s = reverse(s[end-n_species+1:end])
    s_vec = [parse(Int64, ss) for ss in split(s, "")]

    return s_vec
end

function state_idx_to_str(n_species::Int64, idx::Int64)
    #= Get state string from index (single temperature).
        
    Args:
        n_species: Number of species in the system
        idx: Raw index of state
    
    Return:
        str: Raw string of state
    =#
    # Get vectorized state
    species_from_idx = findall(x -> x>0, state_idx_to_vec(n_species, idx))
    str = join(species_from_idx, "*")

    return str
end

function state_idx_to_str(n_species::Int64, 
    n_temps::Int64, 
    idx::Int64)
    #= Get state string from index (multi temperature).
        
    Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
        idx: Raw index of state
    
    Return:
        str: Raw string of state
    =#
    # Get the indices
    (idx_species, idx_temp) = split_idx(n_species, n_temps, idx)

    # Get string from sub function and append temperature
    str = state_idx_to_str(n_species, idx_species) * "|" * string(idx_temp)

    return str
end

function _number_to_base(n::Int64, 
    base::Int64, 
    digits::Int64)
    #= Convert decimal to desired base
        
    Args:
        n: Decimal number
        base: Desired base
        digits: Number of digits in the base
    
    Return:
        n_base: Number in new base
    =#
    n_base = Int64[]
    while n > 0
        push!(n_base, Int64(n % base))
        n ÷= base
    end
    
    while length(n_base) < digits
        push!(n_base, 0)
    end

    return n_base
end

function _action_idx_to_vec(n_species::Int64,
    idx::Int64)
    #= Convert action index to vector (single temperature).
        
    Args:
        n_species: Number of species in the system
        idx: Raw index of action
    
    Return:
        str: Vector of action
    =#
    if idx >= 3^n_species
        throw(DomainError(idx, "This action state index is too large."))
    end
    return _number_to_base(idx, 3, n_species) .- 1
end

function action_idx_to_lv_action(n_species::Int64, 
    n_temps::Int64, 
    idx::Int64)
    #= Convert action index to vector (multi temperature).
        
    Args:
        n_species: Number of species in the system
        n_temps: Number of temperature settings in the system
        idx: Raw index of action
    
    Return:
        lv_action: LVAction object
    =#
    # Get the state and temperature index for the action
    ds_idx = (idx - 1) % (3^n_species)
    dt = ((idx - 1) ÷ (3^n_species)) + 1
    ds = _action_idx_to_vec(n_species, ds_idx)

    # Throw exception
    if dt > n_temps
        throw(DomainError(idx, "This action temperature index is too large."))
    end

    # Additions
    adds = findall(x -> x>0, ds)
    add_string = ""
    if length(adds) > 0
        add_string = "(+" * string(adds)[2:end-1] * ")"
        add_string = replace(add_string, "," => " ")
    end

    # Deletions
    dels = findall(x -> x<0, ds)
    del_string = ""
    if length(dels) > 0
        del_string = "(-" * string(dels)[2:end-1] * ")"
        del_string = replace(del_string, "," => " ")
    end

    # Process and input into dictionary
    # We don't add temperature strings since that is context dependent
    str = add_string * del_string  

    return LVAction(ds, dt, str)
end