using XLSX

# Column order matching DataWorkflow.jl / load_input_data expectations
const _SCENARIO_COLS = [
    "ID", "RCP", "folder", "Purpose", "Intervention", "date", "Region",
    "Year_Start", "Year_End", "init_cover", "fts", "HeatToleranceGroups",
    "HeatToleranceInit", "Heritability", "Plasticity", "counterfactual",
    "Disturbance_file", "Connectivity_file", "Spatial_file", "Geometry_file",
    "Growth_Surv_file", "Fogging_reducer", "Fogging_sites", "Fogging_start",
    "TradeOff", "Deployment area", "TotalCorals", "species", "Enhancement",
    "InterventionYears_start", "duration", "frequency", "Reef_siteids",
    "temp_growth", "species_proportions"
]

# =============================================================================
# FILE DISCOVERY
# =============================================================================

"""
    _discover_data_files(data_dir; spatial_pattern, disturbance_pattern, connectivity_pattern)

Scan `data_dir` and return all valid `(spatial, disturbance, connectivity, geometry)`
file-name quadruples. Files are matched by shared prefix (e.g. `"reef_sites_1"`) so
each spatial is only paired with its own disturbance/connectivity variants.
"""
function _discover_data_files(data_dir::String;
    spatial_pattern::Regex     = r"^reef_sites_.*_nogeo\.RData$",
    disturbance_pattern::Regex = r"temporal_simulation.*\.RData$",
    connectivity_pattern::Regex = r"connectivity.*\.RData$"
)
    files = readdir(data_dir)
    spatial_files = filter(f -> occursin(spatial_pattern, f), files)
    isempty(spatial_files) && error("No spatial files matching $spatial_pattern in $data_dir")

    T = NamedTuple{(:spatial, :disturbance, :connectivity, :geometry), NTuple{4,String}}
    triples = T[]
    for sf in spatial_files
        prefix   = replace(sf, "_nogeo.RData" => "")
        dist_fns = filter(f -> startswith(f, prefix) && occursin(disturbance_pattern, f), files)
        conn_fns = filter(f -> startswith(f, prefix) && occursin(connectivity_pattern, f), files)
        isempty(dist_fns) && @warn "No disturbance files for prefix '$prefix'"
        isempty(conn_fns) && @warn "No connectivity files for prefix '$prefix'"
        for d in dist_fns, c in conn_fns
            push!(triples, (spatial=sf, disturbance=d, connectivity=c, geometry=prefix*".RData"))
        end
    end
    @info "$(length(triples)) file triple(s) from $(length(spatial_files)) spatial file(s)"
    return triples
end

# =============================================================================
# XLSX I/O
# =============================================================================

function _write_scenario_table(df::DataFrame, fpath::String)
    filepath = joinpath(fpath, "ScenarioID.xlsx")
    XLSX.writetable(filepath, df; overwrite=true, sheetname="ScenarioID")
    @info "Written $(nrow(df)) scenario rows to $filepath"
end

function _read_scenario_table(fpath::String)::DataFrame
    filepath = joinpath(fpath, "ScenarioID.xlsx")
    return DataFrame(XLSX.readtable(filepath, "ScenarioID"; infer_eltypes=true))
end

# =============================================================================
# DATAFRAME BUILDERS
# =============================================================================

function _counterfactual_dataframe(file_triples;
    rcps, plasticity_values, year_start_values, year_end_values,
    init_cover_values, heat_tolerance_init_values, heritability_values,
    fts, growth_surv_files, heat_tolerance_groups, region, folder, purpose, temp_growth
)
    today = Dates.format(Dates.today(), "yyyy/mm/dd")

    # Cartesian product — tuple order: (triple, rcp, plasticity, yr_start, yr_end,
    #                                    init_cov, ht_init, heritability)
    combos = vec(collect(Iterators.product(
        file_triples, rcps, plasticity_values,
        year_start_values, year_end_values, init_cover_values,
        heat_tolerance_init_values, heritability_values
    )))
    n = length(combos)

    df = DataFrame(
        ID                      = 1:n,
        RCP                     = [c[2] for c in combos],
        folder                  = fill(folder, n),
        Purpose                 = fill(purpose, n),
        Intervention            = fill("No", n),
        date                    = fill(today, n),
        Region                  = fill(region, n),
        Year_Start              = [c[4] for c in combos],
        Year_End                = [c[5] for c in combos],
        init_cover              = [c[6] for c in combos],
        fts                     = fill(fts, n),
        HeatToleranceGroups     = fill(heat_tolerance_groups, n),
        HeatToleranceInit       = [c[7] for c in combos],
        Heritability            = [c[8] for c in combos],
        Plasticity              = [c[3] for c in combos],
        counterfactual          = fill(missing, n),
        Disturbance_file        = [c[1].disturbance for c in combos],
        Connectivity_file       = [c[1].connectivity for c in combos],
        Spatial_file            = [c[1].spatial for c in combos],
        Geometry_file           = [c[1].geometry for c in combos],
        Growth_Surv_file        = fill(growth_surv_files, n),
        Fogging_reducer         = fill(missing, n),
        Fogging_sites           = fill(missing, n),
        Fogging_start           = fill(missing, n),
        TradeOff                = fill(missing, n),
        TotalCorals             = fill(missing, n),
        species                 = fill(missing, n),
        Enhancement             = fill(missing, n),
        InterventionYears_start = fill(missing, n),
        duration                = fill(missing, n),
        frequency               = fill(missing, n),
        Reef_siteids            = fill(missing, n),
        temp_growth             = fill(temp_growth, n),
        species_proportions     = fill(missing, n)
    )
    df[!, "Deployment area"] = fill(missing, n)
    select!(df, _SCENARIO_COLS)
    return df
end

function _intervention_dataframe(ranking::DataFrame, file_triples;
    n_sites_options, site_selections,
    rcps, year_start_values, year_end_values, init_cover_values,
    heat_tolerance_init_values, heritability_values,
    deployment_areas, total_corals_multipliers, species,
    species_proportions_values, enhancements, intervention_years_start,
    durations, frequencies,
    fts, growth_surv_files, heat_tolerance_groups, region, folder, purpose,
    temp_growth, plasticity, id_start::Int
)
    today   = Dates.format(Dates.today(), "yyyy/mm/dd")
    n_total = nrow(ranking)

    unknown = setdiff(site_selections, ("top", "bottom"))
    isempty(unknown) || @warn "Unknown site selections skipped: $unknown"
    valid_selections = filter(s -> s in ("top", "bottom"), site_selections)

    # Build (n_sites, selection, reef_siteids) configs
    site_configs = vec(map(
        collect(Iterators.product(n_sites_options, valid_selections))
    ) do (n_sites, sel)
        ids = sel == "top" ?
            ranking.site_id[ranking.rank .<= n_sites] :
            ranking.site_id[ranking.rank .> (n_total - n_sites)]
        (n_sites=n_sites, selection=sel, reef_siteids=join(ids, "/"))
    end)

    # Cartesian product — tuple order:
    # (site_cfg, triple, rcp, yr_start, yr_end, init_cov, ht_init, heritability,
    #  dep_area, mult, sp_prop, enhancement, int_year, duration, frequency)
    combos = vec(collect(Iterators.product(
        site_configs, file_triples, rcps,
        year_start_values, year_end_values, init_cover_values,
        heat_tolerance_init_values, heritability_values,
        deployment_areas, total_corals_multipliers, species_proportions_values,
        enhancements, intervention_years_start, durations, frequencies
    )))
    n = length(combos)

    df = DataFrame(
        ID                      = id_start:(id_start + n - 1),
        RCP                     = [c[3]  for c in combos],
        folder                  = fill(folder, n),
        Purpose                 = fill(purpose, n),
        Intervention            = fill("Yes", n),
        date                    = fill(today, n),
        Region                  = fill(region, n),
        Year_Start              = [c[4]  for c in combos],
        Year_End                = [c[5]  for c in combos],
        init_cover              = [c[6]  for c in combos],
        fts                     = fill(fts, n),
        HeatToleranceGroups     = fill(heat_tolerance_groups, n),
        HeatToleranceInit       = [c[7]  for c in combos],
        Heritability            = [c[8]  for c in combos],
        Plasticity              = fill(plasticity, n),
        counterfactual          = fill(missing, n),
        Disturbance_file        = [c[2].disturbance for c in combos],
        Connectivity_file       = [c[2].connectivity for c in combos],
        Spatial_file            = [c[2].spatial for c in combos],
        Geometry_file           = [c[2].geometry for c in combos],
        Growth_Surv_file        = fill(growth_surv_files, n),
        Fogging_reducer         = fill(missing, n),
        Fogging_sites           = fill(missing, n),
        Fogging_start           = fill(missing, n),
        TradeOff                = fill(missing, n),
        TotalCorals             = [c[10] * c[9] for c in combos],
        species                 = fill(species, n),
        Enhancement             = [c[12] for c in combos],
        InterventionYears_start = [c[13] for c in combos],
        duration                = [c[14] for c in combos],
        frequency               = [c[15] for c in combos],
        Reef_siteids            = [c[1].reef_siteids for c in combos],
        temp_growth             = fill(temp_growth, n),
        species_proportions     = [c[11] for c in combos]
    )
    df[!, "Deployment area"] = [c[9] for c in combos]
    select!(df, _SCENARIO_COLS)
    return df
end

# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================

"""
    setup_counterfactuals(fpath, fun_path; kwargs...) -> Union{CScapeResultSet, DataFrame}

**Step 1** of the two-step scenario table workflow.

Auto-discovers spatial × disturbance × connectivity file combinations in `fpath/data/`,
cross-products them with any additional parameter vectors, writes counterfactual rows to
`ScenarioID.xlsx`, runs the simulations, and returns the loaded `CScapeResultSet`.

Set `run=false` to only write the table without running.

# Required keyword arguments
- `fts`: Functional type names, slash-separated (e.g. `"acro_table/acro_corym/..."`)
- `growth_surv_files`: Demography file names per FT, slash-separated

# Sweep parameters (each Vector → one row per value in the cross-product)
- `rcps`, `plasticity_values`, `year_start_values`, `year_end_values`
- `init_cover_values`, `heat_tolerance_init_values`, `heritability_values`
"""
function setup_counterfactuals(fpath::String, fun_path::String;
    fts::String,
    growth_surv_files::String,
    heat_tolerance_groups::String       = "-5_8_14",
    region::String                      = "",
    folder::String                      = "scenarios",
    purpose::String                     = "",
    temp_growth::Int                    = 1,
    rcps::Vector{Int}                   = [1],
    plasticity_values::Vector{Float64}  = [0.0],
    year_start_values::Vector{Int}      = [2025],
    year_end_values::Vector{Int}        = [2050],
    init_cover_values::Vector{String}   = ["0.1_0.1_0.1_0.1_0.1_0.1"],
    heat_tolerance_init_values::Vector{String} = ["0_1.91"],
    heritability_values::Vector{String} = ["0.5_0.01"],
    spatial_pattern::Regex              = r"^reef_sites_.*_nogeo\.RData$",
    disturbance_pattern::Regex          = r"temporal_simulation.*\.RData$",
    connectivity_pattern::Regex         = r"connectivity.*\.RData$",
    n_workers::Int                      = max(1, Sys.CPU_THREADS - 1),
    run::Bool                           = true
)
    file_triples = _discover_data_files(joinpath(fpath, "data");
        spatial_pattern, disturbance_pattern, connectivity_pattern)

    df = _counterfactual_dataframe(file_triples;
        rcps, plasticity_values, year_start_values, year_end_values,
        init_cover_values, heat_tolerance_init_values, heritability_values,
        fts, growth_surv_files, heat_tolerance_groups, region, folder, purpose, temp_growth)

    _write_scenario_table(df, fpath)
    @info "Counterfactual table: $(nrow(df)) rows"

    run || return df

    run_cscape_parallel(collect(df.ID), fpath, fun_path;
        n_workers, export_adria=false, calc_indicators=true)

    return load_results(CScapeResultSet, fpath)
end


"""
    setup_interventions(fpath; ranking, n_counterfactuals, kwargs...) -> DataFrame

**Step 2** of the two-step scenario table workflow.

Takes a pre-computed MCDA `ranking` DataFrame (from `compute_rankings`), selects
top/bottom N sites, cross-products with all intervention parameters, and appends
intervention rows to the existing `ScenarioID.xlsx`.

# Required keyword arguments
- `ranking`: DataFrame with columns `site_id` and `rank` (output of `compute_rankings`)
- `n_counterfactuals`: ID of the last counterfactual row; interventions start at `n+1`
- `fts`, `growth_surv_files`: same values used in `setup_counterfactuals`

# Intervention sweep parameters
- `n_sites_options`, `site_selections` (`"top"` / `"bottom"`)
- `deployment_areas`, `total_corals_multipliers`, `species_proportions_values`
- `enhancements`, `intervention_years_start`, `durations`, `frequencies`
"""
function setup_interventions(fpath::String;
    ranking::DataFrame,
    n_counterfactuals::Int,
    fts::String,
    growth_surv_files::String,
    heat_tolerance_groups::String       = "-5_8_14",
    region::String                      = "",
    folder::String                      = "scenarios",
    purpose::String                     = "",
    temp_growth::Int                    = 1,
    plasticity::Float64                 = 0.0,
    rcps::Vector{Int}                   = [1],
    year_start_values::Vector{Int}      = [2025],
    year_end_values::Vector{Int}        = [2050],
    init_cover_values::Vector{String}   = ["0.1_0.1_0.1_0.1_0.1_0.1"],
    heat_tolerance_init_values::Vector{String} = ["0_1.91"],
    heritability_values::Vector{String} = ["0.5_0.01"],
    spatial_pattern::Regex              = r"^reef_sites_.*_nogeo\.RData$",
    disturbance_pattern::Regex          = r"temporal_simulation.*\.RData$",
    connectivity_pattern::Regex         = r"connectivity.*\.RData$",
    file_triples::Union{Nothing,Vector} = nothing,
    n_sites_options::Vector{Int}        = [5],
    site_selections::Vector{String}     = ["top", "bottom"],
    deployment_areas::Vector{Float64}   = [5000.0],
    total_corals_multipliers::Vector{Float64} = [1.0],
    species::String                     = "1_2_3_4_5",
    species_proportions_values::Vector{String} = ["0.2_0.2_0.2_0.2_0.2"],
    enhancements::Vector{String}        = ["3_0.87"],
    intervention_years_start::Vector{Int} = [2026],
    durations::Vector{Int}              = [5],
    frequencies::Vector{Int}            = [5]
)
    file_triples = if isnothing(file_triples)
        _discover_data_files(joinpath(fpath, "data");
            spatial_pattern, disturbance_pattern, connectivity_pattern)
    else
        file_triples
    end

    int_df = _intervention_dataframe(ranking, file_triples;
        n_sites_options, site_selections,
        rcps, year_start_values, year_end_values, init_cover_values,
        heat_tolerance_init_values, heritability_values,
        deployment_areas, total_corals_multipliers, species,
        species_proportions_values, enhancements, intervention_years_start,
        durations, frequencies,
        fts, growth_surv_files, heat_tolerance_groups, region, folder, purpose,
        temp_growth, plasticity, id_start=n_counterfactuals + 1)

    full_df = _read_scenario_table(fpath)
    allowmissing!(full_df)
    append!(full_df, int_df)

    _write_scenario_table(full_df, fpath)
    @info "$(nrow(int_df)) intervention rows appended ($(nrow(full_df)) total)"

    return int_df
end


"""
    build_scenario_table(fpath, fun_path; kwargs...) -> DataFrame

Combined workflow: runs `setup_counterfactuals`, computes MCDA rankings from the
first counterfactual result, then calls `setup_interventions`.

Returns the complete scenario `DataFrame` (counterfactuals + interventions).

Accepts all keyword arguments from both `setup_counterfactuals` and
`setup_interventions`. Use `setup_counterfactuals` / `setup_interventions`
separately when you need to inspect results or supply a custom ranking between steps.

# Example
```julia
setup_r_environment(fun_path)
df = build_scenario_table(fpath, fun_path;
    fts             = "acro_table/acro_corym/corym_non_acro/small_massive/large_massive/brooder",
    growth_surv_files = "demog_inputs_acro_table.rds/.../demog_inputs_brooder.rds",
    rcps            = [1, 2],
    year_start_values = [2025, 2030],
    deployment_areas  = [5000.0, 10000.0],
    n_sites_options   = [5, 10],
    site_selections   = ["top", "bottom"]
)
```
"""
function build_scenario_table(fpath::String, fun_path::String;
    fts::String,
    growth_surv_files::String,
    heat_tolerance_groups::String       = "-5_8_14",
    region::String                      = "",
    folder::String                      = "scenarios",
    purpose::String                     = "",
    temp_growth::Int                    = 1,
    rcps::Vector{Int}                   = [1],
    plasticity_values::Vector{Float64}  = [0.0],
    plasticity::Float64                 = 0.0,
    year_start_values::Vector{Int}      = [2025],
    year_end_values::Vector{Int}        = [2050],
    init_cover_values::Vector{String}   = ["0.1_0.1_0.1_0.1_0.1_0.1"],
    heat_tolerance_init_values::Vector{String} = ["0_1.91"],
    heritability_values::Vector{String} = ["0.5_0.01"],
    spatial_pattern::Regex              = r"^reef_sites_.*_nogeo\.RData$",
    disturbance_pattern::Regex          = r"temporal_simulation.*\.RData$",
    connectivity_pattern::Regex         = r"connectivity.*\.RData$",
    n_sites_options::Vector{Int}        = [5],
    site_selections::Vector{String}     = ["top", "bottom"],
    deployment_areas::Vector{Float64}   = [5000.0],
    total_corals_multipliers::Vector{Float64} = [1.0],
    species::String                     = "1_2_3_4_5",
    species_proportions_values::Vector{String} = ["0.2_0.2_0.2_0.2_0.2"],
    enhancements::Vector{String}        = ["3_0.87"],
    intervention_years_start::Vector{Int} = [2026],
    durations::Vector{Int}              = [5],
    frequencies::Vector{Int}            = [5],
    n_workers::Int                      = max(1, Sys.CPU_THREADS - 1),
    mcda_prefs::Dict                    = MCDA_PREFS
)
    shared = (fts=fts, growth_surv_files=growth_surv_files,
              heat_tolerance_groups=heat_tolerance_groups,
              region=region, folder=folder, purpose=purpose, temp_growth=temp_growth)
    # year_end_values is kept separate so counterfactuals and interventions can differ
    swept_no_end = (rcps=rcps, year_start_values=year_start_values,
                    init_cover_values=init_cover_values,
                    heat_tolerance_init_values=heat_tolerance_init_values,
                    heritability_values=heritability_values)
    pats   = (spatial_pattern=spatial_pattern, disturbance_pattern=disturbance_pattern,
              connectivity_pattern=connectivity_pattern)

    # Counterfactuals only need to run up to the year before the first intervention
    cf_year_end = [minimum(intervention_years_start) - 1]

    setup_counterfactuals(fpath, fun_path;
        shared..., swept_no_end..., year_end_values=cf_year_end, pats...,
        plasticity_values, n_workers)

    n_counter = nrow(_read_scenario_table(fpath))
    output    = load_output(fpath, 1)
    ranking   = compute_rankings(output, mcda_prefs)

    setup_interventions(fpath;
        ranking, n_counterfactuals=n_counter,
        shared..., swept_no_end..., year_end_values=year_end_values, pats...,
        plasticity, n_sites_options, site_selections,
        deployment_areas, total_corals_multipliers, species,
        species_proportions_values, enhancements, intervention_years_start,
        durations, frequencies)

    return _read_scenario_table(fpath)
end
