"""
    CscapeInterface

Julia wrapper for C-scape coral reef simulation.

# Setup
```julia
using Pkg
Pkg.activate("path/to/C_scape/julia")
Pkg.add("ModelParameters")
Pkg.develop(path="path/to/ADRIAIndicators.jl")
Pkg.instantiate()
```

# Quick Start (Dict - simple)
```julia
include("src/CscapeInterface.jl")
using .CscapeInterface

setup_r_environment()
run_cscape(1, "/path/to/data")
```

# With ModelParameters - Using Helper Functions
```julia
model = load_params(1, fpath)

# Modify parameters individually
set_rcp!(model, 2)
set_year_end!(model, 2050)
set_plasticity!(model, 0.5)
set_dhw_enhance!(model, 3.0)

run_cscape(model)
```

# With ModelParameters - Without Helper Functions
```julia
model = load_params(1, fpath)

# View current values
model[:val]  # (scenario_id, Cyclone_rep, rcp, year_start, year_end, Plasticity, DHW_enhance, draw)

# Method 1: Set all values at once
model[:val] = (1, 1, 2, 2024, 2050, 0.5, 3.0, 0)

# Method 2: Modify single value using collect/Tuple
v = collect(model[:val])
v[3] = 2      # rcp at index 3
v[5] = 2050   # year_end at index 5
model[:val] = Tuple(v)

# Method 3: Using splatting
v = model[:val]
model[:val] = (v[1:2]..., 2, v[4:8]...)  # change rcp (index 3) to 2

run_cscape(model)
```

# Parameter Index Reference
| Index | Parameter   | Description           |
|-------|-------------|-----------------------|
| 1     | scenario_id | Scenario ID           |
| 2     | Cyclone_rep | Cyclone replicate     |
| 3     | rcp         | RCP scenario (1-4)    |
| 4     | year_start  | Start year            |
| 5     | year_end    | End year              |
| 6     | Plasticity  | Plasticity (0-1)      |
| 7     | DHW_enhance | DHW enhancement       |
| 8     | draw        | Posterior draw (0=mean)|
"""
module CscapeInterface

using RCall
using DataFrames
using Statistics
using ModelParameters
using ModelParameters: Model, Param, stripparams, parent
using Distributed
using CSV
using Printf
using TimerOutputs
using Dates
using YAXArrays
using DimensionalData: Dim
using JLD2
using WGLMakie
import ADRIA


# Include source files
include("params.jl")           # ModelParameters types
include("r_wrapper.jl")        # Call R, modify inputs
include("data_access.jl")      # Load outputs, yearly access
include("indicators.jl")       # ADRIAIndicators integration
include("cscape_result_io.jl") # Multi-scenario ResultSet
include("dynamic_ranking.jl")  # Dynamic MCDA re-ranking workflow
include("scenario_setup.jl")   # ScenarioID table construction
include("adria_analysis.jl")   # PAWN + ADRIA sensitivity workflow

# Exports - Types
export CscapeParams, CscapeOutput, CScapeResultSet

# Exports - ModelParameters
export Model, Param  # Re-export from ModelParameters
export load_params, to_dict, from_dict
export set_param!, set_scenario_id!, set_cyclone_rep!, set_rcp!
export set_year_start!, set_year_end!, set_plasticity!, set_dhw_enhance!, set_draw!
using ModelParameters: update!
export update!

# Exports - String/Bool/Vector parameter helpers
export set_region!, set_simulation_name!, set_rootdir_data!
export set_spatial_file!, set_connectivity_file!, set_disturbance_file!, set_growth_surv_file!
export set_init_cover!, set_heat_tolerance!, set_heat_init!, set_fts!
export set_tradeoff!, set_output!, set_use_cached_ipm!, set_temp_growth_switch!
export set_intervention!, set_intervened_sites!, set_sites!
export set_params!

       
# Exports - Dict interface 
export setup_r_environment
export load_input_data, modify_input_data!, load_reef_spatial, modify_reef_spatial!
export run_cscape, run_cscape_modified, setup_interventions
export process_cscape_outputs, run_cscape_parallel
export setup_counterfactuals, build_scenario_table

# Exports - Per-timestep control
export initialise_simulation, run_single_year!, run_years!
export get_simulation_state, print_simulation_state, get_site_cover
export get_simulation_output, get_raw_simulation_output, get_combined_cover_for_ranking
export modify_simulation_state!, finalise_simulation

# Exports - Output access
export load_output, load_ranking_inputs, build_cscape_output, get_yearly_data, get_site_data, get_cover_timeseries
export yearly_iterator, site_iterator

# Exports - Indicators
export to_adria_format, calculate_indicators, indicator_summary, export_for_adria, save_indicators

# Exports - ResultSet
export load_results, load_grouped_results, list_groups, save_results
export scenario_outcome, scenario_groups

# Exports - Dynamic reranking workflow
export run_dynamic_reranking   # was main()
export run_shared_phase
export run_mcda_phase
export build_iteration_configs
export define_worker_code!
export compute_rankings
export compute_rankings_from_cover
export MCDA_PREFS              # if kept as module-level const


# Exports - ADRIA sensitivity workflow
export run_adria_pawn

# Re-export ADRIA ResultSet interface helpers (not exported by ADRIA itself)
import ADRIA: n_locations, n_scenarios, loc_k_area, loc_k
export n_locations, n_scenarios, loc_k_area, loc_k

end