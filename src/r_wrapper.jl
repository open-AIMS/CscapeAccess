# =============================================================================
# r_wrapper.jl - Julia-R Interface for C-scape
# =============================================================================
#
# This module provides Julia functions to run C-scape simulations via RCall.
# Julia handles: IPM batch creation, simulation orchestration, intervention setup
#
# Main entry points:
#   - run_cscape(scenario_id, fpath)     → Run full simulation
#   - initialise_simulation(...)          → Year-by-year control
#   - run_single_year!(env, year)         → Run one year
#   - modify_simulation_state!(env; ...)  → Modify parameters mid-run
#
# =============================================================================


# =============================================================================
# SETUP FUNCTIONS
# =============================================================================

"""
    setup_r_environment(fun_path::String; enable_parallel::Bool=true, n_cores::Int=0)

Setup R environment by sourcing all C-scape functions.

# Arguments
- `fun_path::String`: Path to C_scape directory (default: parent of julia/)
- `enable_parallel::Bool`: Enable parallel processing (default: true)
- `n_cores::Int`: Number of cores for parallel (default: 0 = auto-detect)

# Example
```julia
setup_r_environment()  # Auto-detect path, parallel enabled
setup_r_environment(enable_parallel=false)  # No parallel
setup_r_environment(n_cores=4)  # Use 4 cores
```
"""
function setup_r_environment(fun_path::String = ""; enable_parallel::Bool = true, n_cores::Int = 0)
    # Auto-detect path if not provided
    if isempty(fun_path)
        fun_path = dirname(dirname(@__DIR__))  # Go up from julia/src/ to C_scape/
    end
    
    @info "Setting up R environment from: $fun_path"
    
    R"""
    # Free large globals from prior runs in the persistent RCall session.
    stale_objects <- c("MainEnvir", "year_output", "sim_output", "output",
                       "BATCH_IPM_CACHE", "BATCH_IPM_N_SITES", "create_ipm_original")
    existing_objects <- intersect(stale_objects, ls(envir = .GlobalEnv))
    if (length(existing_objects) > 0) {
        rm(list = existing_objects, envir = .GlobalEnv)
    }
    invisible(gc())

    setwd($fun_path)
    
    # Tell main.R we're calling from Julia:
    #   - Skips R version check in setup_cscape()
    #   - Skips library(JuliaCall) to avoid circular dependency
    FROM_JULIA <- TRUE
    source('main.R')
    
    # Source C-scape function files
    source('cscape_sim.R')
    source('intervention_setup.R')
    source('ancillary_functions.R')
    source('annual_site_taxa_loops.R')
    source('ipm_pred.R')
    
    # Source modular functions if available
    modular_path <- file.path($fun_path, "modules")
    if (dir.exists(modular_path)) {
        modular_files <- list.files(modular_path, pattern = "\\.R$", full.names = TRUE)
        for (f in modular_files) source(f)
        cat("Loaded", length(modular_files), "modular functions\n")
    }
    
    # Source all toolbox functions
    toolbox_files <- list.files(file.path($fun_path, "1_toolbox"), pattern = "\\.R$", full.names = TRUE)
    for (f in toolbox_files) {
      source(f)
    }
    
    # Store path for later use
    .fun_path <- $fun_path
    
    print(paste("R environment ready!", length(toolbox_files), "toolbox functions loaded"))
    """
    
    # Enable parallel if requested
    if enable_parallel
        @rput n_cores
        R"""
        ENABLE_PARALLEL <- TRUE
        library(doParallel)
        library(foreach)
        if ($n_cores == 0) {
            n_cores_use <- parallel::detectCores() - 1
        } else {
            n_cores_use <- $n_cores
        }
        registerDoParallel(cores = n_cores_use)
        print(paste("Parallel enabled:", n_cores_use, "cores"))
        """
        @info "Parallel processing enabled"
    else
        R"""
        ENABLE_PARALLEL <- FALSE
        """
        @info "Parallel processing disabled"
    end
    
    @info "R environment setup complete"
    return fun_path
end


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

"""
    load_input_data(scenario_id::Int, fpath::String) -> Dict

Load InputData from ScenarioID.xlsx as Julia Dict.

# Example
```julia
input_data = load_input_data(1, "/path/to/data")
input_data["year_end"] = 2050  # Modify in Julia
run_cscape(input_data)
```
"""
function load_input_data(scenario_id::Int, fpath::String)
    
    @rput scenario_id fpath
    
    R"""
    # Create log file
    log_dir <- paste0($fpath, "/logs")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- paste0(log_dir, "/juvenile_cap_log", $scenario_id, ".csv")
    
    # Load scenarios
    Scenarios <- read_excel(paste0($fpath, "/ScenarioID.xlsx"), 
                            sheet="ScenarioID", na = c("NA", "N/A", "na"))
    scn_row <- Scenarios[Scenarios$ID == $scenario_id, ]
    
    # Parse fields
    fts <- strsplit(as.character(scn_row$fts), split="/")[[1]]
    demog <- strsplit(as.character(scn_row$Growth_Surv_file), split="/")[[1]]
    
    # Build InputData
    InputData_R <- list(
        rootdir_data = $fpath,
        scenario_id = $scenario_id,
        region = scn_row$Region,
        simulation_name = scn_row$folder,
        draw = NA,
        Cyclone_rep = 1,
        rcp = scn_row$RCP,
        Intervention = scn_row$Intervention,
        sites = TRUE,
        year_start = as.integer(scn_row$Year_Start),
        year_end = as.integer(scn_row$Year_End),
        init_cover = scn_row$init_cover,
        fts = fts,
        HeatTolerance = scn_row$HeatToleranceGroups,
        HeatInit = scn_row$HeatToleranceInit,
        Heritability = scn_row$Heritability,
        Plasticity = scn_row$Plasticity,
        output = TRUE,
        TradeOff = scn_row$TradeOff,
        DHW_enhance = as.numeric(scn_row$Enhancement),
        Spatial_file = scn_row$Spatial_file,
        Connectivity_file = scn_row$Connectivity_file,
        Disturbance_file = scn_row$Disturbance_file,
        Growth_Surv_file = demog,
        intervened_sites = scn_row$Reef_siteids,
        log_file = log_file,
        use_cached_IPM = TRUE,
        temp_growth_switch = scn_row$temp_growth
    )
    """
    
    input_data_raw = rcopy(R"InputData_R")
    
    # Convert to Dict{String, Any} for easier access
    input_data = Dict{String, Any}(string(k) => v for (k, v) in pairs(input_data_raw))
    
    @info "Loaded InputData for scenario $scenario_id" region=input_data["region"] years="$(input_data["year_start"])-$(input_data["year_end"])"
    
    return input_data
end


"""
    modify_input_data!(input_data::Dict; kwargs...)

Modify InputData in place.

# Example
```julia
input_data = load_input_data(1, fpath)
modify_input_data!(input_data; year_end=2050, draw=5, Plasticity=0.6)
run_cscape(input_data)
```
"""
function modify_input_data!(input_data::Dict; kwargs...)
    for (key, value) in kwargs
        key_str = string(key)
        old_val = get(input_data, key_str, nothing)
        input_data[key_str] = value
        @info "Modified $key_str: $old_val → $value"
    end
    return input_data
end


"""
    setup_interventions(scenario_id::Int, fpath::String)

Setup intervention files required before running simulation.
Called automatically by run_cscape() and initialise_simulation().
"""
function setup_interventions(scenario_id::Int, fpath::String)
    
    @rput scenario_id fpath
    
    R"""
    Scenarios <- read_excel(paste0($fpath, "/ScenarioID.xlsx"), 
                            sheet="ScenarioID", na = c("NA", "N/A", "na"))
    scn <- $scenario_id
    
    demog <- strsplit(as.character(Scenarios$Growth_Surv_file[Scenarios$ID==scn]), split="/")[[1]]
    
    # Coral deployment
    Deployment_Sites <- strsplit(as.character(Scenarios$Reef_siteids[Scenarios$ID==scn]), split="/")[[1]]
    Deployment_area <- Scenarios$`Deployment area`[Scenarios$ID==scn]
    Year_start <- Scenarios$InterventionYears_start[Scenarios$ID==scn]
    Duration <- Scenarios$duration[Scenarios$ID==scn]
    Frequency <- Scenarios$frequency[Scenarios$ID==scn]
    Enhancement <- Scenarios$Enhancement[Scenarios$ID==scn]
    Species <- as.numeric(strsplit(as.character(Scenarios$species[Scenarios$ID==scn]), split="_")[[1]])
    Proportions <- as.numeric(strsplit(as.character(Scenarios$species_proportions[Scenarios$ID == scn]), split = "_")[[1]])
  
    TotalCorals <- Scenarios$TotalCorals[Scenarios$ID==scn]
    
    if (!is.na(Year_start)) {
    if (!is.na(Deployment_area)){
          CoralAddition <- setup_coral_addition(Deployment_Sites, Deployment_area, TotalCorals, Year_start,
                                          Duration, Frequency, Enhancement, Species,Proportions, scn, fpath, demog)
    } else {
          CoralAddition <- readRDS(paste(fpath, "/data/CoralAddition",scn,".RDS", sep = ""))
    } 
} else {
        CoralAddition <- data.frame(scenario=numeric(), reef_siteid=factor(),
                                    Year=integer(), ft=numeric(),
                    no_int_corals=numeric(), proportion=numeric(), m2=numeric(), density=numeric(),
                                    meshpt_int_corals=character(), Enhancement=numeric())
    }
    
    # Fogging
    Reducer <- Scenarios$Fogging_reducer[Scenarios$ID==scn]
    Sites <- Scenarios$Fogging_sites[Scenarios$ID==scn]
    Start <- Scenarios$Fogging_start[Scenarios$ID==scn]
    
    if (!is.na(Sites)) {
        Fogging <- setup_cooling_shading(Reducer, Sites, Start, scn)
    } else {
        Fogging <- data.frame(scenario=numeric(), reef_siteid=character(),
                              year=numeric(), reduction=numeric())
    }
    
    Interventions <- list(Coral=CoralAddition, Fogging=Fogging)
    
    # Save to temp file first to avoid OneDrive sync issues
    temp_file <- tempfile(fileext = ".RData")
    dir.create(dirname(temp_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(Interventions, file = temp_file)
    
    # Copy to final location
    final_file <- paste0($fpath, "/data/Interventions", scn, ".RData")
    file.copy(temp_file, final_file, overwrite = TRUE)
    file.remove(temp_file)
    
    print(paste("Interventions saved for scenario", scn))
    """
    
    @info "Interventions setup complete for scenario $scenario_id"
end


# =============================================================================
# INDICATOR CALCULATION HELPER
# =============================================================================

"""
    _calculate_and_save_indicators(fpath, scenario_id, draw_val; output, filepath) -> Union{CscapeIndicators, Nothing}

Internal: calculate ADRIA indicators from saved output and save as JLD2.
Called automatically by `run_cscape` and `finalise_simulation`.

`filepath` overrides the default save path. When `finalise_simulation` is called with
a custom `filename`, the indicator file is placed alongside the output file so that
parallel workers running the same scenario_id/draw do not overwrite each other.
"""
function _calculate_and_save_indicators(fpath::String, scenario_id::Int, draw_val;
                                         output = nothing,
                                         filepath::Union{String,Nothing} = nothing)
    draw_str = if isnothing(draw_val) || ismissing(draw_val)
        "NA"
    elseif draw_val isa Number && isnan(draw_val)
        "NA"
    else
        string(Int(draw_val))
    end

    try
        if isnothing(output)
            output = load_output(fpath, scenario_id; draw=draw_str)
        end
        indicators = calculate_indicators(output)

        indicator_path = if isnothing(filepath)
            joinpath(fpath, "adria_exports",
                     "Indicators_scenario_$(scenario_id)_draw_$(draw_str).jld2")
        else
            filepath
        end
        @info "Saving indicators to $indicator_path"
        save_indicators(indicators, output, indicator_path)
        indicator_summary(indicators)
        return indicators
    catch e
        @warn "Could not calculate ADRIA indicators" exception=(e, catch_backtrace())
        return nothing
    end
end


"""
    process_cscape_outputs(fpath; scenario_ids, export_adria, calc_indicators, overwrite) -> CScapeResultSet

Process pre-existing CScape out_array RDS files in `fpath/model_outputs/` without
re-running the simulation. Generates indicator and/or ADRIA export JLD2 files for each
discovered file, then loads and returns a `CScapeResultSet`.

Scenarios whose indicator file already exists in `adria_exports/` are skipped unless
`overwrite=true`.

# Arguments
- `fpath`: Data directory with standard structure (`model_outputs/`, `adria_exports/`,
  `data/`, `ScenarioID.xlsx`)
- `scenario_ids`: Optional filter; default processes all RDS files found
- `export_adria`: Save full `CscapeOutput` to `adria_*.jld2` (large files; default `false`)
- `calc_indicators`: Compute and save `Indicators_*.jld2` (default `true`)
- `overwrite`: Recompute even if the indicator file already exists (default `false`)

# Example
```julia
setup_r_environment(fun_path)
rs = process_cscape_outputs("/path/to/data")
rs = process_cscape_outputs("/path/to/data"; scenario_ids=1:10, export_adria=true)
```
"""
function process_cscape_outputs(fpath::String;
    scenario_ids::Union{Nothing,AbstractVector{Int}} = nothing,
    export_adria::Bool = false,
    calc_indicators::Bool = true,
    overwrite::Bool = false
)
    output_dir = joinpath(fpath, "model_outputs")
    adria_dir  = joinpath(fpath, "adria_exports")
    mkpath(adria_dir)

    rds_pattern = r"Array_scenario_(\d+)_draw_(.+)\.rds"
    all_files = readdir(output_dir)
    matched = filter(f -> occursin(rds_pattern, f), all_files)

    scenarios = map(matched) do f
        m = match(rds_pattern, f)
        (scenario_id=parse(Int, m[1]), draw=String(m[2]))
    end

    if !isnothing(scenario_ids)
        scenarios = filter(s -> s.scenario_id ∈ scenario_ids, scenarios)
    end

    if !overwrite
        scenarios = filter(scenarios) do s
            ind_file = joinpath(adria_dir, "Indicators_scenario_$(s.scenario_id)_draw_$(s.draw).jld2")
            !isfile(ind_file)
        end
    end

    @info "Processing $(length(scenarios)) scenario(s) from $output_dir"

    for s in scenarios
        try
            output = load_output(fpath, s.scenario_id; draw=s.draw)

            if calc_indicators
                indicators = calculate_indicators(output)
                ind_path = joinpath(adria_dir, "Indicators_scenario_$(s.scenario_id)_draw_$(s.draw).jld2")
                @info "Saving indicators to $ind_path"
                save_indicators(indicators, output, ind_path)
                indicator_summary(indicators)
            end

            if export_adria
                adria_path = joinpath(adria_dir, "adria_scenario_$(s.scenario_id)_draw_$(s.draw).jld2")
                @info "Exporting ADRIA output to $adria_path"
                export_for_adria(output, adria_path)
            end
        catch e
            @warn "Failed to process scenario $(s.scenario_id) draw $(s.draw)" exception=(e, catch_backtrace())
        end
    end

    return load_results(CScapeResultSet, fpath)
end


"""
    run_cscape_parallel(scenario_ids, fpath, fun_path; n_workers, export_adria, calc_indicators) -> Vector{Int}

Run multiple C-scape scenarios in parallel using `pmap`. Spawns Julia workers, loads
`CscapeInterface` on each, runs the scenarios, then cleans up workers.

# Arguments
- `scenario_ids`: IDs to run (e.g. `1:100` or `[2, 5, 7]`)
- `fpath`: Data directory passed to `run_cscape`
- `fun_path`: C_scape source directory passed to `setup_r_environment`
- `n_workers`: Number of worker processes (default: `Sys.CPU_THREADS - 1`)
- `export_adria`: Passed to `run_cscape` (default `true`)
- `calc_indicators`: Passed to `run_cscape` (default `true`)

# Returns
Vector of scenario IDs that failed (empty if all succeeded).

# Example
```julia
failed = run_cscape_parallel(1:100, fpath, fun_path)
```
"""
function run_cscape_parallel(scenario_ids, fpath::String, fun_path::String;
    n_workers::Int = max(1, Sys.CPU_THREADS - 1),
    export_adria::Bool = true,
    calc_indicators::Bool = true
)
    addprocs(n_workers)
    try
        @sync for p in workers()
            @async remotecall_wait(Core.eval, p, Main, :(using CscapeInterface))
        end

        @info "Running $(length(scenario_ids)) scenario(s) across $n_workers worker(s)"

        failed = pmap(scenario_ids) do sid
            result = try
                setup_r_environment(fun_path; enable_parallel=false)
                run_cscape(sid, fpath; export_adria=export_adria, calc_indicators=calc_indicators)
                nothing
            catch e
                @error "Scenario $sid failed" exception=e
                sid
            end
            result
        end |> x -> filter(!isnothing, x)

        isempty(failed) ? @info("All scenarios complete") : @warn("Failed scenarios: $failed")
        return failed
    finally
        rmprocs(workers())
    end
end


# =============================================================================
# MAIN SIMULATION FUNCTIONS
# =============================================================================

"""
    run_cscape(fpath::String; kwargs...)
    run_cscape(scenario_id::Int, fpath::String; kwargs...)
    run_cscape(scenario_ids::Vector{Int}, fpath::String; kwargs...)
    run_cscape(input_data::Dict; kwargs...)

Run C-scape simulation.

# Arguments
- `fpath::String`: Data path — runs ALL scenarios from ScenarioID.xlsx
- `scenario_id::Int`: Single scenario ID
- `scenario_ids::Vector{Int}`: Multiple scenario IDs to run sequentially
- `input_data::Dict`: Pre-built InputData dictionary. If `scenario_id` is a vector, 
   runs each ID using the same custom settings.

# Keywords
- `fun_path::String`: Path to C_scape (auto-detected if empty)
- `export_adria::Bool`: Export for ADRIAIndicators after completion (default: true)
- `calc_indicators::Bool`: Calculate and save ADRIA indicators (default: true)

# Examples
```julia
# Run ALL scenarios in ScenarioID.xlsx
run_cscape(fpath)

# Single scenario
run_cscape(1, fpath)

# Multiple scenarios by ID
run_cscape([1, 2, 3, 4, 5, 56], fpath)

# Customised single scenario
input_data = load_input_data(1, fpath)
input_data["year_end"] = 2050
run_cscape(input_data)

# Customised with multiple scenario IDs
input_data = load_input_data(1, fpath)
input_data["scenario_id"] = [1, 2, 3, 4, 5, 56]
input_data["year_end"] = 2050
input_data["Plasticity"] = 0.6
run_cscape(input_data)

# Skip indicator calculation
run_cscape(1, fpath; calc_indicators=false)
```
"""
function run_cscape(fpath::String;
                    fun_path::String = "",
                    export_adria::Bool = true,
                    calc_indicators::Bool = true)
    
    # Read all scenario IDs from ScenarioID.xlsx
    @rput fpath
    R"""
    library(readxl)
    Scenarios <- read_excel(paste0($fpath, "/ScenarioID.xlsx"), 
                            sheet="ScenarioID", na = c("NA", "N/A", "na"))
    all_ids <- as.integer(Scenarios$ID)
    """
    all_ids = Int.(rcopy(R"all_ids"))
    
    @info "Found $(length(all_ids)) scenarios in ScenarioID.xlsx: $all_ids"
    
    return run_cscape(all_ids, fpath; fun_path = fun_path, 
                      export_adria = export_adria, 
                      calc_indicators = calc_indicators)
end


function run_cscape(scenario_id::Int, fpath::String; 
                    fun_path::String = "", 
                    export_adria::Bool = true,
                    calc_indicators::Bool = true)
    
    # Load InputData and run
    input_data = load_input_data(scenario_id, fpath)
    _run_cscape_single(input_data; fun_path = fun_path, 
                       export_adria = export_adria, 
                       calc_indicators = calc_indicators)
end


function run_cscape(scenario_ids::Vector{Int}, fpath::String;
                    fun_path::String = "",
                    export_adria::Bool = true,
                    calc_indicators::Bool = true)
    
    n = length(scenario_ids)
    @info "Running $n scenarios: $scenario_ids"
    
    failed = Int[]
    
    for (i, sid) in enumerate(scenario_ids)
        @info "━━━ Scenario $sid ($i/$n) ━━━"
        try
            run_cscape(sid, fpath; fun_path = fun_path, 
                       export_adria = export_adria, 
                       calc_indicators = calc_indicators)
        catch e
            @error "Scenario $sid failed" exception=(e, catch_backtrace())
            push!(failed, sid)
        end
    end
    
    if !isempty(failed)
        @warn "$(length(failed))/$n scenarios failed: $failed"
    else
        @info "All $n scenarios completed successfully"
    end
    
    return failed
end


function run_cscape(input_data::Dict; 
                    fun_path::String = "",
                    export_adria::Bool = true,
                    calc_indicators::Bool = true)
    
    sid = input_data["scenario_id"]
    
    # If scenario_id is a vector, run each with the same custom settings
    if sid isa AbstractVector
        scenario_ids = Int.(sid)
        fpath = input_data["rootdir_data"]
        n = length(scenario_ids)
        @info "Running $n scenarios with custom settings: $scenario_ids"
        
        failed = Int[]
        
        for (i, id) in enumerate(scenario_ids)
            @info "━━━ Scenario $id ($i/$n) ━━━"
            try
                # Load fresh scenario data, then apply custom overrides
                id_data = load_input_data(id, fpath)
                for (k, v) in input_data
                    if k != "scenario_id"
                        id_data[k] = v
                    end
                end
                id_data["scenario_id"] = id
                
                _run_cscape_single(id_data; fun_path = fun_path,
                                   export_adria = export_adria,
                                   calc_indicators = calc_indicators)
            catch e
                @error "Scenario $id failed" exception=(e, catch_backtrace())
                push!(failed, id)
            end
        end
        
        if !isempty(failed)
            @warn "$(length(failed))/$n scenarios failed: $failed"
        else
            @info "All $n scenarios completed successfully"
        end
        
        return failed
    else
        # Single scenario
        _run_cscape_single(input_data; fun_path = fun_path,
                           export_adria = export_adria,
                           calc_indicators = calc_indicators)
    end
end


# --- Internal: runs exactly one scenario from a fully-prepared Dict ---
function _run_cscape_single(input_data::Dict; 
                            fun_path::String = "",
                            export_adria::Bool = true,
                            calc_indicators::Bool = true)
    
    # Get fun_path from R (set by setup_r_environment)
    if isempty(fun_path)
        fun_path = rcopy(R".fun_path")
    end
    
    scenario_id = input_data["scenario_id"]
    fpath = input_data["rootdir_data"]
    
    @info "Running C-scape" scenario=scenario_id years="$(input_data["year_start"])-$(input_data["year_end"])"
    
    # Setup interventions
    setup_interventions(scenario_id, fpath)
    
    # Pass to R
    @rput input_data fun_path
    
    R"""
    # Source the Julia batch runner
    source(file.path(.fun_path, "run_with_julia_batch.R"))

    # Convert Julia Dict to R list
    InputData <- as.list(input_data)

    # Handle NA
    if (is.null(InputData$draw)) InputData$draw <- NA

    # Run simulation with Julia batch acceleration
    print("Starting simulation...")
    tic()
    MainEnvir <- run_simulation_with_julia_batch(InputData, $fun_path, RubbleHandle = FALSE)
    toc()
    print("Simulation complete!")
    """
    
    # Load output once and share between ADRIA export and indicator calculation
    if export_adria || calc_indicators
        draw_val = get(input_data, "draw", nothing)
        draw_str = if isnothing(draw_val) || ismissing(draw_val) || (draw_val isa Number && isnan(draw_val))
            "NA"
        else
            string(Int(draw_val))
        end

        built_output = try
            load_output(fpath, scenario_id; draw=draw_str, from_envir=true)
        catch e
            @warn "Could not load simulation output from R environment" exception=(e, catch_backtrace())
            nothing
        end

        if export_adria && !isnothing(built_output)
            try
                adria_path = joinpath(fpath, "adria_exports",
                                      "adria_scenario_$(scenario_id)_draw_$(draw_str).jld2")
                export_for_adria(built_output, adria_path)
            catch e
                @warn "ADRIA export failed" exception=(e, catch_backtrace())
            end
        end

        if calc_indicators
            _calculate_and_save_indicators(fpath, scenario_id, draw_val; output=built_output)
        end
    end
    
    @info "Simulation finished!"
end


# =============================================================================
# YEAR-BY-YEAR CONTROL FUNCTIONS
# =============================================================================

"""
    initialise_simulation(scenario_id::Int, fpath::String; fun_path::String="")
    initialise_simulation(input_data::Dict; fun_path::String="")

Initialise simulation and return environment handle for year-by-year control.

# Example
```julia
# Initialise
env = initialise_simulation(1, fpath)

# Run year by year
for year in 2008:2018
    output = run_single_year!(env, year)
    println("Year \$year complete")
end

# Modify parameters mid-run
modify_simulation_state!(env; kappa_scale=0.9, plasticity=0.6)

# Continue running
run_years!(env, 2019:2025)

# Save results
finalise_simulation(env)
```
"""
function initialise_simulation(scenario_id::Int, fpath::String; fun_path::String = "")
    input_data = load_input_data(scenario_id, fpath)
    return initialise_simulation(input_data; fun_path = fun_path)
end

function initialise_simulation(input_data::Dict; fun_path::String = "")
    # Auto-detect fun_path
    if isempty(fun_path)
        fun_path = dirname(dirname(@__DIR__))
    end
    
    # Setup interventions first
    scenario_id = input_data["scenario_id"]
    fpath = input_data["rootdir_data"]
    setup_interventions(scenario_id, fpath)
    
    @rput input_data fun_path
    
    R"""
    # Ensure repeated initialisation does not retain the previous simulation state.
    stale_objects <- c("MainEnvir", "year_output", "sim_output", "output",
                       "BATCH_IPM_CACHE", "BATCH_IPM_N_SITES", "create_ipm_original")
    existing_objects <- intersect(stale_objects, ls(envir = .GlobalEnv))
    if (length(existing_objects) > 0) {
        rm(list = existing_objects, envir = .GlobalEnv)
    }
    invisible(gc())

    # Source initialise function (it sources all other required files)
    source(file.path($fun_path, "modules/run_year.R"))
    
    # Convert Julia Dict to R list
    InputData <- as.list(input_data)
    if (is.null(InputData$draw)) InputData$draw <- NA
    
    # Initialise simulation with Julia batch enabled
    MainEnvir <- initialise_simulation(InputData, $fun_path, RubbleHandle = FALSE, use_julia_batch = TRUE)
    
    # Store fun_path in R global for later use
    .fun_path <- $fun_path
    """
    
    @info "Simulation initialised" scenario=scenario_id years="$(input_data["year_start"])-$(input_data["year_end"])"
    
    # Return a handle (the actual MainEnvir stays in R)
    return Dict(
        "scenario_id" => scenario_id,
        "fpath" => fpath,
        "fun_path" => fun_path,
        "year_start" => input_data["year_start"],
        "year_end" => input_data["year_end"],
        "current_year" => input_data["year_start"],
        "initialised" => true
    )
end


"""
    run_single_year!(env::Dict, year::Int; save_yearly::Bool=false)

Run a single year of simulation.

# Arguments
- `env`: Environment handle from `initialise_simulation`
- `year`: Year to run
- `save_yearly`: Save output for this year to file

# Returns
- Array slice for this year [sites, intervention, ft, enhancement, metrics]

# Example
```julia
env = initialise_simulation(1, fpath)

for year in 2008:2018
    output = run_single_year!(env, year)
    mean_cover = mean(output[:, 1, :, :, 1])
    println("Year \$year: mean cover = \$mean_cover")
end
```
"""
function run_single_year!(env::Dict, year::Int; save_yearly::Bool = false)
    if !get(env, "initialised", false)
        error("Simulation not initialised. Call initialise_simulation() first.")
    end
    
    fun_path = env["fun_path"]
    @rput year save_yearly fun_path
    
    R"""
    # run_single_year uses Julia batch internally if enabled during init
    year_output <- run_single_year($year, MainEnvir, save_yearly = $save_yearly)
    """
    
    year_output = rcopy(R"year_output")
    env["current_year"] = year
    
    @info "Completed year $year"
    
    return year_output
end


"""
    run_years!(env::Dict, years::UnitRange{Int}; save_yearly::Bool=false)

Run multiple years of simulation.

# Example
```julia
env = initialise_simulation(1, fpath)
run_years!(env, 2008:2012)  # Run first 5 years

# Analyze or modify state here...

run_years!(env, 2013:2018)  # Run remaining years
```
"""
function run_years!(env::Dict, years::UnitRange{Int}; save_yearly::Bool = false)
    outputs = []
    for year in years
        output = run_single_year!(env, year; save_yearly = save_yearly)
        push!(outputs, output)
    end
    return outputs
end


"""
    get_simulation_output(env::Dict)

Get the full output array from simulation.

# Returns
- Full output array [years, sites, intervention, ft, enhancement, metrics]
"""
function get_simulation_output(env::Dict)
    R"""
    out_array_full <- MainEnvir$out_array
    """
    return rcopy(R"out_array_full")
end


"""
    get_raw_simulation_output(env::Dict) -> (Array, Vector{Int}, Vector{String}, Vector{String})

Extract raw simulation array and dimension names directly from R memory (no disk I/O).
Returns `(raw_array, years, site_ids, fts)`.

`MainEnvir\$out_array` always has 2 intervention dimensions; the combined slot is added
by `build_cscape_output` in Julia. Aliases `raw_arr` / `dimnames_arr` are freed
immediately after transfer — `MainEnvir\$out_array` is unaffected (R copy-on-modify).
"""
function get_raw_simulation_output(env::Dict)
    R"""
    raw_arr      <- MainEnvir$out_array
    dimnames_arr <- dimnames(raw_arr)
    """
    raw_array  = rcopy(R"raw_arr")
    dimnames_r = rcopy(R"dimnames_arr")
    R"rm(raw_arr, dimnames_arr); invisible(gc())"

    years    = parse.(Int, dimnames_r[:year])
    site_ids = String.(dimnames_r[:reef_siteid])
    fts      = String.(dimnames_r[:ft])
    return raw_array, years, site_ids, fts
end


"""
    get_combined_cover_for_ranking(env::Dict, area_mat::Matrix{Float64}) -> Array{Float64,4}

Extract the cover metric (metric index 1) from `MainEnvir\$out_array` and return the
area-weighted combined cover as `[years, sites, n_ft, n_enh]`.

Transfers only ~1/106 of the full array from R; no CscapeOutput is constructed.
`MainEnvir\$out_array` is unaffected (R copy-on-modify on the cover_slice alias).
"""
function get_combined_cover_for_ranking(env::Dict, area_mat::Matrix{Float64})
    R"""
    cover_slice <- MainEnvir$out_array[, , , , , 1]
    """
    cover_raw = rcopy(R"cover_slice")     # [years, sites, 2, n_ft, n_enh]
    R"rm(cover_slice); invisible(gc())"

    n_years, n_sites, _, n_ft, n_enh = size(cover_raw)
    combined = Array{Float64,4}(undef, n_years, n_sites, n_ft, n_enh)
    for site in 1:n_sites
        tot      = area_mat[site, 3]
        prop_non = tot > 0.0 ? area_mat[site, 1] / tot : 1.0
        prop_int = tot > 0.0 ? area_mat[site, 2] / tot : 0.0
        @views combined[:, site, :, :] .=
            prop_non .* cover_raw[:, site, 1, :, :] .+
            prop_int .* cover_raw[:, site, 2, :, :]
    end
    return combined
end


"""
    get_simulation_state(env::Dict)

Get current simulation state information including all modifiable parameters.

# Returns
Dict with current_year, n_sites, n_fts, kappa stats, plasticity, heritability, etc.
"""
function get_simulation_state(env::Dict)
    R"""
    suppressWarnings({
    # Helper to safely get stats
    safe_mean <- function(x) {
      if (is.null(x) || length(x) == 0) return(NA)
      if (is.numeric(x)) return(mean(x, na.rm = TRUE))
      return(NA)
    }
    safe_min <- function(x) {
      if (is.null(x) || length(x) == 0) return(NA)
      if (is.numeric(x)) return(min(x, na.rm = TRUE))
      return(NA)
    }
    safe_max <- function(x) {
      if (is.null(x) || length(x) == 0) return(NA)
      if (is.numeric(x)) return(max(x, na.rm = TRUE))
      return(NA)
    }
    
    # Get connectivity stats
    connec_mean <- NA
    connec_max <- NA
    tryCatch({
      if (!is.null(MainEnvir$connec)) {
        if (is.matrix(MainEnvir$connec)) {
          connec_mean <- mean(MainEnvir$connec, na.rm = TRUE)
          connec_max <- max(MainEnvir$connec, na.rm = TRUE)
        } else if (is.list(MainEnvir$connec)) {
          numeric_vals <- c()
          for (item in MainEnvir$connec) {
            if (is.numeric(item)) {
              numeric_vals <- c(numeric_vals, as.vector(item))
            } else if (is.matrix(item)) {
              numeric_vals <- c(numeric_vals, as.vector(item))
            }
          }
          if (length(numeric_vals) > 0) {
            connec_mean <- mean(numeric_vals, na.rm = TRUE)
            connec_max <- max(numeric_vals, na.rm = TRUE)
          }
        }
      }
    }, error = function(e) {
      connec_mean <- NA
      connec_max <- NA
    })
    
    # Get dhw_threshold
    dhw_thresh <- MainEnvir$dhw_threshold
    if (is.null(dhw_thresh) || !is.numeric(dhw_thresh)) {
      dhw_thresh <- MainEnvir$reef_spatial$bleach_threshold
    }
    
    sim_state <- list(
        current_year = MainEnvir$current_year,
        n_sites = length(MainEnvir$site_names),
        n_fts = length(MainEnvir$fts),
        fts = MainEnvir$fts,
        site_names_sample = MainEnvir$site_names[1:min(5, length(MainEnvir$site_names))],
        years = c(MainEnvir$InputData$year_start, MainEnvir$InputData$year_end),
        
        # Kappa stats
        kappa_mean = safe_mean(MainEnvir$kappa),
        kappa_min = safe_min(MainEnvir$kappa),
        kappa_max = safe_max(MainEnvir$kappa),
        
        # DHW threshold stats
        dhw_threshold_mean = safe_mean(dhw_thresh),
        dhw_threshold_min = safe_min(dhw_thresh),
        dhw_threshold_max = safe_max(dhw_thresh),
        
        # Modifiable parameters
        plasticity = MainEnvir$Plasticity,
        heritability = MainEnvir$Heritability,
        dhw_enhance = MainEnvir$InputData$DHW_enhance,
        
        # Reef and connectivity
        reef_areas_total = sum(MainEnvir$reef_areas, na.rm = TRUE),
        reef_areas_mean = mean(MainEnvir$reef_areas, na.rm = TRUE),
        connectivity_mean = connec_mean,
        connectivity_max = connec_max,
        
        # Interventions
        rubble_handle = MainEnvir$RubbleHandle,
        fogging_sites = if (!is.null(MainEnvir$Fogging)) nrow(MainEnvir$Fogging) else 0,
        intervention_sites = if (!is.null(MainEnvir$CoralIntervention)) length(unique(MainEnvir$CoralIntervention$reef_siteid)) else 0
    )
    })
    """
    return rcopy(R"sim_state")
end


"""
    print_simulation_state(env::Dict)

Print a formatted summary of current simulation state.
"""
function print_simulation_state(env::Dict)
    state = get_simulation_state(env)
    
    # Helper to safely format numbers
    fmt(x, digits=3) = ismissing(x) || isnothing(x) || (x isa Number && isnan(x)) ? "NA" : round(x, digits=digits)
    
    # Helper to get value with either Symbol or String key
    getval(d, k) = haskey(d, k) ? d[k] : (haskey(d, Symbol(k)) ? d[Symbol(k)] : "NA")
    
    println("\n" * "="^50)
    println("SIMULATION STATE")
    println("="^50)
    println("Year: $(getval(state, "current_year")) | Range: $(getval(state, "years")[1])-$(getval(state, "years")[2])")
    println("Sites: $(getval(state, "n_sites")) | FTs: $(getval(state, "n_fts"))")
    println("-"^50)
    println("MODIFIABLE PARAMETERS:")
    println("  Kappa:         mean=$(fmt(getval(state, "kappa_mean"))), min=$(fmt(getval(state, "kappa_min"))), max=$(fmt(getval(state, "kappa_max")))")
    println("  DHW threshold: mean=$(fmt(getval(state, "dhw_threshold_mean"))), min=$(fmt(getval(state, "dhw_threshold_min"))), max=$(fmt(getval(state, "dhw_threshold_max")))")
    println("  Plasticity:    $(getval(state, "plasticity"))")
    println("  Heritability:  $(getval(state, "heritability"))")
    println("  DHW enhance:   $(getval(state, "dhw_enhance"))")
    println("  Reef areas:    total=$(fmt(getval(state, "reef_areas_total"), 2)), mean=$(fmt(getval(state, "reef_areas_mean"), 4))")
    println("  Connectivity:  mean=$(fmt(getval(state, "connectivity_mean"), 6)), max=$(fmt(getval(state, "connectivity_max"), 6))")
    println("  Rubble handle: $(getval(state, "rubble_handle"))")
    println("-"^50)
    println("INTERVENTIONS:")
    println("  Fogging sites:      $(getval(state, "fogging_sites"))")
    println("  Intervention sites: $(getval(state, "intervention_sites"))")
    println("="^50 * "\n")
    
    return state
end


"""
    get_site_cover(env::Dict) -> Dict

Get current coral cover and carrying capacity for all sites.
Use this to identify sites with available space before deploying corals.

# Returns
Dict with:
- `site_names`: Vector of site IDs
- `cover`: Total coral cover per site
- `kappa`: Carrying capacity per site
- `remaining_capacity`: kappa - cover per site
- `proportion_full`: cover / kappa per site (0-1)

# Example
```julia
env = initialise_simulation(1, fpath)
run_years!(env, 2008:2012)

# Check site availability
state = get_site_cover(env)

# Find sites with >20% remaining capacity
available = findall(i -> state["proportion_full"][i] < 0.8, 
                    1:length(state["site_names"]))
println("Sites with space: ", state["site_names"][available])

# Deploy corals only at available sites
deploy_sites = state["site_names"][available[1:min(5, length(available))]]
```
"""
function get_site_cover(env::Dict)
    R"""
    n_sites <- length(MainEnvir$site_names)
    n_fts <- length(MainEnvir$fts)
    year_idx <- MainEnvir$current_year - MainEnvir$InputData$year_start + 1
    
    # Calculate total cover per site from out_array
    # out_array dimensions: [year, site, intervention, ft, enhancement, metrics]
    # metric index 1 = coral cover
    site_cover <- rep(0, n_sites)
    
    for (s in 1:n_sites) {
        total <- 0
        for (ft in 1:n_fts) {
            # Sum cover across interventions and enhancements
            cover_vals <- MainEnvir$out_array[year_idx, s, , ft, , 1]
            total <- total + sum(cover_vals, na.rm = TRUE)
        }
        site_cover[s] <- total
    }
    
    kappa_vals <- MainEnvir$kappa
    remaining <- kappa_vals - site_cover
    remaining[remaining < 0] <- 0
    
    prop_full <- ifelse(kappa_vals > 0, site_cover / kappa_vals, 1)
    prop_full[prop_full > 1] <- 1
    
    cover_state <- list(
        site_names = MainEnvir$site_names,
        cover = site_cover,
        kappa = kappa_vals,
        remaining_capacity = remaining,
        proportion_full = prop_full
    )
    """
    return rcopy(R"cover_state")
end


"""
    modify_simulation_state!(env::Dict; kwargs...)

Modify simulation state between years.

# Keywords
- `kappa_scale::Float64`: Scale all kappa values (e.g., 0.9 = reduce by 10%)
- `kappa_values::Vector{Float64}`: Set kappa values directly
- `dhw_threshold::Vector{Float64}`: Set DHW threshold directly per site
- `plasticity::Float64`: Set plasticity value (0-1)
- `heritability::String`: Set heritability as "mean_sd" (e.g., "0.3_0.01")
- `dhw_enhance::Float64`: Set DHW enhancement value
- `connectivity_scale::Float64`: Scale connectivity matrix values
- `fogging_reduction::Float64`: Set fogging DHW reduction factor
- `rubble_handle::Bool`: Enable/disable rubble handling
- `coral_deployment::Dict`: Replace coral deployment schedule (DataFrame as Dict)
    Required keys: reef_siteid, Year, ft, no_int_corals, proportion, m2, density, meshpt_int_corals, Enhancement
- `fogging_schedule::Dict`: Replace fogging schedule (DataFrame as Dict)
  Required keys: reef_siteid, year, reduction
- `cots_mortality::Dict`: Modify COTS values in reef_temporal
  Keys: year (Int), site_indices (Vector{Int}), values (Matrix 4×n_sites for bins 1-4)
  Or: year (Int), scale (Float64) to scale all COTS for that year

# Example
```julia
env = initialise_simulation(1, fpath)
run_years!(env, 2008:2012)

# Check which sites have space for coral deployment
cover_state = get_site_cover(env)
available = findall(i -> cover_state["remaining_capacity"][i] > 1000, 1:length(cover_state["remaining_capacity"]))

# Modify deployment schedule based on site availability
modify_simulation_state!(env; 
    kappa_scale=0.9,
    coral_deployment=Dict(
        "reef_siteid" => ["Moore_001", "Moore_002"],
        "Year"        => [2013, 2013],
        "ft"          => [1, 1],
        "no_int_corals" => [1000.0, 1000.0],
        "proportion"  => [1.0, 1.0],
        "m2"          => [500.0, 500.0],
        "density"     => [2.0, 2.0],
        "meshpt_int_corals" => ["4", "4"],
        "Enhancement" => [3, 3]
    ),
    fogging_schedule=Dict(
        "reef_siteid" => ["Moore_001", "Moore_002"],
        "year"        => [2013, 2013],
        "reduction"   => [0.3, 0.3]
    ),
    cots_mortality=Dict(
        "year" => 2013,
        "scale" => 0.5    # halve COTS pressure for year 2013
    ) # or
    cots_mortality=Dict(
        "year"         => 2014,
        "site_indices" => [1],
        "values"       => [0.5; 0.3; 0.1; 0.05]  # 4 bins for 1 site (column vector)
    ) 
)

run_years!(env, 2013:2018)
finalise_simulation(env)
```
"""
function modify_simulation_state!(env::Dict; 
                                   kappa_scale::Union{Float64,Nothing} = nothing,
                                   kappa_values::Union{Vector{Float64},Nothing} = nothing,
                                   dhw_threshold::Union{Vector{Float64},Nothing} = nothing,
                                   plasticity::Union{Float64,String,Nothing} = nothing,
                                   heritability::Union{String,Nothing} = nothing,
                                   dhw_enhance::Union{Float64,Nothing} = nothing,
                                   connectivity_scale::Union{Float64,Nothing} = nothing,
                                   fogging_reduction::Union{Float64,Nothing} = nothing,
                                   rubble_handle::Union{Bool,Nothing} = nothing,
                                   coral_deployment::Union{Dict,Nothing} = nothing,
                                   fogging_schedule::Union{Dict,Nothing} = nothing,
                                   cots_mortality::Union{Dict,Nothing} = nothing)
    
    # Kappa scaling
    if !isnothing(kappa_scale)
        @rput kappa_scale
        R"""
        MainEnvir$kappa <- MainEnvir$kappa * $kappa_scale
        MainEnvir$totalreef_areas <- sum(MainEnvir$reef_areas * MainEnvir$kappa / 100)
        print(paste("Scaled kappa by", $kappa_scale, "| New mean:", round(mean(MainEnvir$kappa), 3)))
        """
        @info "Scaled kappa by $kappa_scale"
    end
    
    # Kappa direct values
    if !isnothing(kappa_values)
        @rput kappa_values
        R"""
        n_sites <- length(MainEnvir$site_names)
        if (length($kappa_values) == n_sites) {
            MainEnvir$kappa <- $kappa_values
            MainEnvir$totalreef_areas <- sum(MainEnvir$reef_areas * MainEnvir$kappa / 100)
            print(paste("Set kappa values | New mean:", round(mean(MainEnvir$kappa), 3)))
        } else {
            stop(paste("kappa_values length", length($kappa_values), "must match number of sites", n_sites))
        }
        """
        @info "Set kappa values directly"
    end
    
    # DHW threshold direct values
    if !isnothing(dhw_threshold)
        @rput dhw_threshold
        R"""
        n_sites <- length(MainEnvir$site_names)
        if (length($dhw_threshold) == n_sites) {
            MainEnvir$dhw_threshold <- $dhw_threshold
            print(paste("Set DHW threshold | New mean:", round(mean(MainEnvir$dhw_threshold, na.rm=TRUE), 3)))
        } else {
            stop(paste("dhw_threshold length", length($dhw_threshold), "must match number of sites", n_sites))
        }
        """
        @info "Set DHW threshold directly"
    end
    
    # Plasticity
    if !isnothing(plasticity)
        @rput plasticity
        R"""
        MainEnvir$Plasticity <- $plasticity
        MainEnvir$InputData$Plasticity <- $plasticity
        print(paste("Set Plasticity to:", $plasticity))
        """
        @info "Set Plasticity to $plasticity"
    end
    
    # Heritability
    if !isnothing(heritability)
        @rput heritability
        R"""
        MainEnvir$Heritability <- $heritability
        MainEnvir$InputData$Heritability <- $heritability
        print(paste("Set Heritability to:", $heritability))
        """
        @info "Set Heritability to $heritability"
    end
    
    # DHW enhance
    if !isnothing(dhw_enhance)
        @rput dhw_enhance
        R"""
        MainEnvir$DHW_enhance <- rep($dhw_enhance, length(MainEnvir$HeatToleranceClasses))
        MainEnvir$InputData$DHW_enhance <- $dhw_enhance
        print(paste("Set DHW_enhance to:", $dhw_enhance))
        """
        @info "Set DHW_enhance to $dhw_enhance"
    end
    
    # Connectivity scaling
    if !isnothing(connectivity_scale)
        @rput connectivity_scale
        R"""
        if (is.list(MainEnvir$connec)) {
          for (i in seq_along(MainEnvir$connec)) {
            if (is.numeric(MainEnvir$connec[[i]])) {
              MainEnvir$connec[[i]] <- MainEnvir$connec[[i]] * $connectivity_scale
            }
          }
          print(paste("Scaled connectivity (list) by", $connectivity_scale))
        } else if (is.matrix(MainEnvir$connec) || is.numeric(MainEnvir$connec)) {
          MainEnvir$connec <- MainEnvir$connec * $connectivity_scale
          print(paste("Scaled connectivity by", $connectivity_scale))
        } else {
          print("Warning: connectivity structure not recognized, not scaled")
        }
        """
        @info "Scaled connectivity by $connectivity_scale"
    end
    
    # Fogging reduction
    if !isnothing(fogging_reduction)
        @rput fogging_reduction
        R"""
        if (!is.null(MainEnvir$Fogging) && nrow(MainEnvir$Fogging) > 0) {
            MainEnvir$Fogging$reduction <- $fogging_reduction
            print(paste("Set Fogging reduction to:", $fogging_reduction))
        } else {
            print("No Fogging data to modify")
        }
        """
        @info "Set Fogging reduction to $fogging_reduction"
    end
    
    # Rubble handle
    if !isnothing(rubble_handle)
        @rput rubble_handle
        R"""
        MainEnvir$RubbleHandle <- $rubble_handle
        print(paste("Set RubbleHandle to:", $rubble_handle))
        """
        @info "Set RubbleHandle to $rubble_handle"
    end
    
    # Coral deployment schedule (replace entire DataFrame)
    if !isnothing(coral_deployment)
        @rput coral_deployment
        R"""
        new_coral <- as.data.frame(coral_deployment, stringsAsFactors = FALSE)
        
        # Validate required columns
        required <- c("reef_siteid", "Year", "ft", "no_int_corals", "proportion", "m2", "density",
                       "meshpt_int_corals", "Enhancement")
        missing_cols <- setdiff(required, names(new_coral))
        if (length(missing_cols) > 0) {
            stop(paste("coral_deployment missing columns:", paste(missing_cols, collapse=", ")))
        }
        
        # Add scenario column if missing
        if (!"scenario" %in% names(new_coral)) {
            new_coral$scenario <- MainEnvir$InputData$scenario_id
        }
        
        MainEnvir$CoralIntervention <- new_coral
        
        # Update intervened sites list
        new_sites <- unique(new_coral$reef_siteid)
        MainEnvir$Intervened_sites <- union(MainEnvir$Intervened_sites, new_sites)
        
        # CRITICAL: Split current population into intervention slot (slot 2)
        # Without this, out_array[..., "yes", ...] is all zeros, causing
        # cots_mortality.R SizeAv to produce vector NaN instead of scalar
        # This mirrors initialise_populations() for intervention scenarios
        current_year <- MainEnvir$current_year
        year_names <- dimnames(MainEnvir$out_array)[[1]]
        year_idx <- which(year_names == as.character(current_year))
        if (length(year_idx) == 0) {
            year_idx <- max(which(as.numeric(year_names) <= current_year))
        }
        
        n_fts <- dim(MainEnvir$out_array)[4]
        reef_areas <- MainEnvir$reef_areas
        
        # Pre-build O(1) lookup maps (replaces O(n_sites) which() per deploy site)
        site_idx_map <- setNames(seq_along(MainEnvir$site_names), MainEnvir$site_names)
        m2_lookup    <- setNames(
            new_coral$m2[match(new_sites, new_coral$reef_siteid)],
            new_sites
        )

        # Unset from environment so arr is the sole owner (NAMED == 1).
        # Assignments to arr[...] then modify in-place -- no full-array copies.
        arr <- MainEnvir$out_array
        MainEnvir$out_array <- NULL

        for (site in new_sites) {
            site_idx <- site_idx_map[[site]]
            if (is.null(site_idx) || is.na(site_idx)) next

            # Skip if intervention slot already has population data
            if (sum(abs(arr[year_idx, site_idx, "yes", , , ]), na.rm = TRUE) > 0) next

            prop <- min(m2_lookup[[site]] / reef_areas[site_idx], 1.0)

            # Extract small per-site slices -- new objects, do not increment NAMED of arr
            no_sz <- arr[year_idx, site_idx, "no", , , 2:106]  # [n_ft, n_enh, 105]
            no_cv <- arr[year_idx, site_idx, "no", , , 1]       # [n_ft, n_enh]

            # Assign in-place (vectorised over all FTs); no FT loop needed
            arr[year_idx, site_idx, "yes", , , 2:106] <- prop       * no_sz
            arr[year_idx, site_idx, "no",  , , 2:106] <- (1 - prop) * no_sz
            arr[year_idx, site_idx, "yes", , , 1]     <- no_cv
            rm(no_sz, no_cv)

            cat("  Split population at", site, "(prop=", round(prop, 3), ")\n")
        }

        # Restore array to environment (one copy)
        MainEnvir$out_array <- arr
        rm(arr)
        
        print(paste("Updated coral deployment:", nrow(new_coral), "entries across", 
                     length(new_sites), "sites"))
        print(paste("  Years:", paste(sort(unique(new_coral$Year)), collapse=", ")))
        print(paste("  Sites:", paste(head(new_sites, 5), collapse=", "), 
                     if(length(new_sites) > 5) "..." else ""))
        """
        @info "Updated coral deployment schedule"
    end
    
    # Fogging schedule (replace entire DataFrame)
    if !isnothing(fogging_schedule)
        @rput fogging_schedule
        R"""
        new_fogging <- as.data.frame(fogging_schedule, stringsAsFactors = FALSE)
        
        # Validate required columns
        required <- c("reef_siteid", "year", "reduction")
        missing_cols <- setdiff(required, names(new_fogging))
        if (length(missing_cols) > 0) {
            stop(paste("fogging_schedule missing columns:", paste(missing_cols, collapse=", ")))
        }
        
        # Add scenario column if missing
        if (!"scenario" %in% names(new_fogging)) {
            new_fogging$scenario <- MainEnvir$InputData$scenario_id
        }
        
        MainEnvir$Fogging <- new_fogging
        
        fogging_sites <- unique(new_fogging$reef_siteid)
        print(paste("Updated fogging schedule:", nrow(new_fogging), "entries across",
                     length(fogging_sites), "sites"))
        print(paste("  Years:", paste(sort(unique(new_fogging$year)), collapse=", ")))
        print(paste("  Reduction range:", 
                     round(min(new_fogging$reduction), 3), "-", 
                     round(max(new_fogging$reduction), 3)))
        """
        @info "Updated fogging schedule"
    end
    
    # COTS mortality modification in reef_temporal
    # disturbances is 3D: [year, site, disturbance_type] where COTS bins are indices 4:7
    if !isnothing(cots_mortality)
        @rput cots_mortality
        R"""
        year_target <- cots_mortality$year
        years_temporal <- as.numeric(dimnames(MainEnvir$reef_temporal$disturbances)$year)
        year_idx <- which(years_temporal == year_target)
        
        if (length(year_idx) == 0) {
            stop(paste("COTS year", year_target, "not found in reef_temporal. Available:", 
                        paste(range(years_temporal), collapse="-")))
        }
        
        if (!is.null(cots_mortality$scale)) {
            # Scale mode: multiply all COTS bins for the target year
            scale_factor <- cots_mortality$scale
            old_mean <- mean(MainEnvir$reef_temporal$disturbances[year_idx, , 4:7], na.rm=TRUE)
            MainEnvir$reef_temporal$disturbances[year_idx, , 4:7] <- 
                MainEnvir$reef_temporal$disturbances[year_idx, , 4:7] * scale_factor
            new_mean <- mean(MainEnvir$reef_temporal$disturbances[year_idx, , 4:7], na.rm=TRUE)
            print(paste("Scaled COTS for year", year_target, "by", scale_factor,
                         "| Mean:", round(old_mean, 4), "->", round(new_mean, 4)))
            
        } else if (!is.null(cots_mortality$site_indices) && !is.null(cots_mortality$values)) {
            # Direct mode: set COTS bins for specific sites
            site_idx <- as.integer(cots_mortality$site_indices)
            vals <- cots_mortality$values  # Expected: 4 x n_sites matrix (bins x sites)
            
            if (!is.matrix(vals)) vals <- matrix(vals, nrow=4)
            if (ncol(vals) != length(site_idx)) {
                stop(paste("COTS values columns", ncol(vals), 
                           "must match site_indices length", length(site_idx)))
            }
            
            for (s in seq_along(site_idx)) {
                MainEnvir$reef_temporal$disturbances[year_idx, site_idx[s], 4:7] <- vals[, s]
            }
            print(paste("Set COTS for year", year_target, "at", length(site_idx), "sites"))
        } else {
            stop("cots_mortality requires either 'scale' or both 'site_indices' and 'values'")
        }
        """
        @info "Modified COTS mortality"
    end
    
    return env
end


"""
    finalise_simulation(env::Dict; export_adria::Bool=true, calc_indicators::Bool=true, filename::Union{String,Nothing}=nothing)

Finalise simulation, save outputs, calculate indicators, and clean up.

# Arguments
- `env`: Simulation environment from initialise_simulation
- `export_adria`: Export for ADRIA format (default: true)
- `calc_indicators`: Calculate and save ADRIA indicators (default: true)
- `filename`: Custom filename or full path. If nothing, uses default naming.

# Examples
```julia
env = initialise_simulation(1, fpath)
for year in 2008:2018
    run_single_year!(env, year)
end

# Default naming (with indicators)
finalise_simulation(env)

# Custom filename
finalise_simulation(env, filename="my_custom_run")

# Skip indicators
finalise_simulation(env, calc_indicators=false)
```
"""
function finalise_simulation(env::Dict; export_adria::Bool = true,
                              calc_indicators::Bool = true,
                              filename::Union{String,Nothing} = nothing,
                              keep_open::Bool = false)
    fpath = env["fpath"]
    scenario_id = env["scenario_id"]
    fun_path = env["fun_path"]

    @rput fpath scenario_id fun_path

    # Capture draw and load output from R env BEFORE calling save_outputs.
    # save_outputs.R may modify MainEnvir$out_array (e.g. adding a combined intervention
    # slot), which would break a subsequent from_envir load.
    draw_val = rcopy(R"MainEnvir$InputData$draw")
    draw_str = if isnothing(draw_val) || ismissing(draw_val) || (draw_val isa Number && isnan(draw_val))
        "NA"
    else
        string(Int(draw_val))
    end

    built_output = if export_adria || calc_indicators
        try
            load_output(fpath, scenario_id; draw=draw_str, from_envir=true)
        catch e
            @warn "Could not load simulation output from R environment" exception=(e, catch_backtrace())
            nothing
        end
    else
        nothing
    end

    if isnothing(filename)
        R"""
        source(file.path($fun_path, "modules/save_outputs.R"))
        save_outputs(MainEnvir$out_array, MainEnvir$InputData, MainEnvir$InputData$draw)
        """
    else
        @rput filename
        R"""
        source(file.path($fun_path, "modules/save_outputs.R"))
        save_outputs(MainEnvir$out_array, MainEnvir$InputData, MainEnvir$InputData$draw, filename = $filename)
        """
    end

    # Derive output-collocated paths when a custom filename is given,
    # so parallel workers with the same scenario_id/draw don't overwrite each other.
    adria_path, indicator_path = if isnothing(filename)
        default_dir = joinpath(fpath, "adria_exports")
        (joinpath(default_dir, "adria_scenario_$(scenario_id)_draw_$(draw_str).jld2"),
         joinpath(default_dir, "Indicators_scenario_$(scenario_id)_draw_$(draw_str).jld2"))
    else
        base_name = splitext(basename(filename))[1]
        out_dir   = dirname(filename) == "." ? joinpath(fpath, "adria_exports") : dirname(filename)
        (joinpath(out_dir, "adria_" * base_name * ".jld2"),
         joinpath(out_dir, "Indicators_" * base_name * ".jld2"))
    end

    if export_adria && !isnothing(built_output)
        @info "Saving ADRIA export to $adria_path"
        try
            export_for_adria(built_output, adria_path)
        catch e
            @warn "ADRIA export failed" exception=(e, catch_backtrace())
        end
    end

    # Calculate and save ADRIA indicators
    if calc_indicators
        _calculate_and_save_indicators(fpath, scenario_id, draw_val;
                                        output=built_output, filepath=indicator_path)
    end
    
    if !keep_open
       env["initialised"] = false
    end
    
    @info "Simulation finalised"
end

"""
    load_params(scenario_id::Int, fpath::String) -> Model{CscapeParams}

Load InputData as ModelParameters Model for ADRIA interoperability.

# Example
```julia
model = load_params(1, fpath)

# View parameters
model           # Shows table
model[:val]     # Get values as tuple
model[:bounds]  # Get bounds

# Modify with helper functions
set_year_end!(model, 2050)
set_plasticity!(model, 0.5)
set_rcp!(model, 2)

# Or modify directly
v = collect(model[:val])
v[5] = 2050  # year_end at index 5
model[:val] = Tuple(v)

# Run
run_cscape(model)
```
"""
function load_params(scenario_id::Int, fpath::String)
    input_data = load_input_data(scenario_id, fpath)
    return from_dict(input_data)
end