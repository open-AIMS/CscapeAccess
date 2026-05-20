"""
    run_adria_pawn(rs, fpath; ts_window=5)

Run the full PAWN sensitivity + ADRIA visualisation workflow on a `CScapeResultSet`.

Outputs are written under `fpath/pawn_inputs/` (CSV + JLD2) and `fpath/figures/` (PNG).

Returns a NamedTuple with fields:
  `input_table`, `pawn_df`, `summary_df`, `y_*` outcome vectors, `pawn_*` sensitivity results.

# Example
```julia
using CscapeInterface
fpath = "/path/to/SingleReef"
rs    = load_results(CScapeResultSet, joinpath(fpath, "All_results.jld2"))
res   = run_adria_pawn(rs, fpath)
```
"""
function run_adria_pawn(
    rs, fpath::AbstractString;
    ts_window::Int = 5,
)
    ts_all = collect(timesteps(rs))

    # ── 1. Split intervention / counterfactual masks ──────────────────────────
    interv_col  = uppercase.(strip.(string.(rs.inputs[!, :Intervention])))
    interv_mask = interv_col .== "YES"
    cf_mask     = interv_col .== "NO"
    interv_idx  = findall(interv_mask)
    cf_idx_all  = findall(cf_mask)

    @info "Scenarios" n_intervention=sum(interv_mask) n_counterfactual=sum(cf_mask)

    # ── 2. Match pairs by shared environmental parameters ─────────────────────
    match_cols = Symbol.(filter(c -> c in names(rs.inputs),
        ["RCP", "Disturbance_file", "Connectivity_file", "Spatial_file"]))

    @info "Matching pairs by: $(match_cols)"

    pair_cf_idx = map(interv_idx) do ii
        row = rs.inputs[ii, :]
        findfirst(cf_idx_all) do ji
            all(string(rs.inputs[ji, col]) == string(row[col]) for col in match_cols)
        end
    end

    unmatched = isnothing.(pair_cf_idx)
    if any(unmatched)
        @warn "$(sum(unmatched)) intervention scenarios had no counterfactual match — dropping"
        interv_idx  = interv_idx[.!unmatched]
        pair_cf_idx = pair_cf_idx[.!unmatched]
    end
    pair_cf_idx = Int.(pair_cf_idx)

    @info "Matched pairs: $(length(interv_idx))"

    # ── 3. Intervention site index helper ────────────────────────────────────
    function interv_site_idx(scenario_row_idx::Int)
        site_id = strip(string(rs.inputs[scenario_row_idx, :Reef_siteids]))
        idx = findfirst(==(site_id), rs.loc_ids)
        isnothing(idx) && error("Reef site '$site_id' not found in rs.loc_ids")
        return idx
    end

    # ── 4. Outcome differences at intervention site ───────────────────────────
    ts_range = max(1, length(ts_all) - ts_window + 1):length(ts_all)

    taxa_arr    = parent(rs.outcomes[:relative_loc_taxa_cover])
    cover_arr   = parent(rs.outcomes[:relative_cover])
    juv_arr     = parent(rs.outcomes[:relative_juveniles])
    fish_arr    = parent(rs.outcomes[:reef_fish_index])
    tourism_arr = parent(rs.outcomes[:reef_tourism_index])
    even_arr    = parent(rs.outcomes[:coral_evenness])
    div_arr     = parent(rs.outcomes[:coral_diversity])
    shelter_arr = parent(rs.outcomes[:relative_shelter_volume])
    biodiv_arr  = parent(rs.outcomes[:reef_biodiversity_condition_index])

    function y_diff_at_site(arr3d, si, ii, cf_ii)
        mean(@view(arr3d[ts_range, si, ii]) .- @view(arr3d[ts_range, si, cf_ii]))
    end

    function y_diff_taxa_at_site(si, ii, cf_ii)
        diff = taxa_arr[ts_range, :, si, ii] .- taxa_arr[ts_range, :, si, cf_ii]
        mean(sum(diff, dims=2))
    end

    y_cover   = Float64[y_diff_at_site(cover_arr,   interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_juv     = Float64[y_diff_at_site(juv_arr,     interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_fish    = Float64[y_diff_at_site(fish_arr,    interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_tourism = Float64[y_diff_at_site(tourism_arr, interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_even    = Float64[y_diff_at_site(even_arr,    interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_div     = Float64[y_diff_at_site(div_arr,     interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_shelter = Float64[y_diff_at_site(shelter_arr, interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_biodiv  = Float64[y_diff_at_site(biodiv_arr,  interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]
    y_taxa    = Float64[y_diff_taxa_at_site(interv_site_idx(ii), ii, pair_cf_idx[k]) for (k,ii) in enumerate(interv_idx)]

    @info "Outcome differences computed" n_pairs=length(y_cover)

    # ── 5. Connectivity matrix ────────────────────────────────────────────────
    spatial   = rs.loc_data
    conn_file = let rel = string(rs.inputs[interv_idx[1], :Connectivity_file])
        f = joinpath(fpath, "data", rel)
        isfile(f) ? f : joinpath(fpath, rel)
    end
    @rput conn_file
    R"""
    conn_raw <- tryCatch(
        readRDS(conn_file),
        error = function(e) {
            env_tmp <- new.env()
            load(conn_file, envir = env_tmp)
            env_tmp[[ls(env_tmp)[1]]]
        }
    )
    conn_mat_r <- if (is.array(conn_raw) && length(dim(conn_raw)) == 3) {
        conn_raw[,,1]
    } else if (is.data.frame(conn_raw)) {
        as.matrix(conn_raw)
    } else {
        as.matrix(conn_raw)
    }
    """
    conn_mat = Matrix{Float64}(rcopy(R"conn_mat_r"))
    @info "Connectivity matrix" size=size(conn_mat)

    # ── 6. Build input table ──────────────────────────────────────────────────
    function disturbance_id(dist_file::String)
        m = match(r"_([^_]+)\.RData$", dist_file)
        isnothing(m) ? dist_file : m.captures[1]
    end

    function parse_enhancement(enh_str::String)
        tokens = split(enh_str, '_')
        isempty(tokens) ? NaN : parse(Float64, tokens[1])
    end

    function parse_proportions(prop_str::String)
        Float64[parse(Float64, t) for t in split(strip(prop_str), '_')]
    end

    dist_id_strings = [disturbance_id(string(rs.inputs[ii, :Disturbance_file])) for ii in interv_idx]
    unique_dist     = unique(dist_id_strings)
    dist_id_map     = Dict(v => Float64(i) for (i, v) in enumerate(unique_dist))

    juv_base_arr  = parent(rs.outcomes[:relative_juveniles])
    even_base_arr = parent(rs.outcomes[:coral_evenness])
    div_base_arr  = parent(rs.outcomes[:coral_diversity])

    input_table = DataFrame()
    for (k, ii) in enumerate(interv_idx)
        si    = interv_site_idx(ii)
        cf_ii = pair_cf_idx[k]
        row   = rs.inputs[ii, :]

        t_start = Int(row[:InterventionYears_start])
        t_pre   = something(findlast(ts_all .< t_start), 1)

        dep_area = Float64(row[Symbol("Deployment area")])
        tot_cor  = Float64(row[:TotalCorals])
        enh      = parse_enhancement(string(row[:Enhancement]))
        sp_props = parse_proportions(string(row[:species_proportions]))

        push!(input_table, Dict(
            :Deployment_area     => dep_area,
            :TotalCorals         => tot_cor,
            :Enhancement         => enh,
            :coral_density       => tot_cor / dep_area,
            :species_proportions => string(row[:species_proportions]),
            :sp_prop_1           => get(sp_props, 1, NaN),
            :sp_prop_2           => get(sp_props, 2, NaN),
            :sp_prop_3           => get(sp_props, 3, NaN),
            :sp_prop_4           => get(sp_props, 4, NaN),
            :sp_prop_5           => get(sp_props, 5, NaN),
            :disturbance_id      => dist_id_strings[k],
            :disturbance_num     => dist_id_map[dist_id_strings[k]],
            :site_k              => spatial[si, :k],
            :site_depth_med      => spatial[si, :depth_med],
            :site_area           => spatial[si, :area],
            :site_hab_area       => spatial[si, :area] * spatial[si, :k],
            :site_ub_med         => spatial[si, :ub_med],
            :conn_self           => conn_mat[si, si],
            :conn_incoming       => sum(conn_mat[:, si]) - conn_mat[si, si],
            :conn_outgoing       => sum(conn_mat[si, :]) - conn_mat[si, si],
            :baseline_cover      => sum(taxa_arr[t_pre, :, si, cf_ii]),
            :baseline_evenness   => even_base_arr[t_pre, si, cf_ii],
            :baseline_diversity  => div_base_arr[t_pre, si, cf_ii],
            :baseline_juveniles  => juv_base_arr[t_pre, si, cf_ii],
        ); cols=:union, promote=true)
    end

    @info "Input table built" rows=nrow(input_table) cols=ncol(input_table)

    # ── 7. Save CSV + JLD2 outputs ────────────────────────────────────────────
    out_dir = joinpath(fpath, "pawn_inputs")
    mkpath(out_dir)
    CSV.write(joinpath(out_dir, "pawn_input_table.csv"), input_table)

    JLD2.jldsave(joinpath(out_dir, "pawn_y_vectors.jld2");
        y_cover, y_juv, y_fish, y_tourism,
        y_even, y_div, y_shelter, y_biodiv, y_taxa,
        interv_scenario_ids = Int.(rs.inputs[interv_idx, :scenario_id]))

    @info "PAWN inputs saved" path=out_dir n_pairs=length(y_cover)

    # ── 8. Summary table ──────────────────────────────────────────────────────
    summary_df = DataFrame(
        scenario_id = Int.(rs.inputs[interv_idx,  :scenario_id]),
        cf_scenario = Int.(rs.inputs[pair_cf_idx, :scenario_id]),
        Δcover      = round.(y_cover,   digits=4),
        Δjuveniles  = round.(y_juv,     digits=4),
        Δfish       = round.(y_fish,    digits=4),
        Δtourism    = round.(y_tourism, digits=4),
        Δtaxa_cover = round.(y_taxa,    digits=4),
    )
    println("\n=== Intervention vs Counterfactual Summary ===")
    println(summary_df)

    # ── 9. PAWN sensitivity analysis ─────────────────────────────────────────
    numeric_cols = filter(c -> eltype(input_table[!, c]) <: Real, names(input_table))
    pawn_df = DataFrames.select(input_table, numeric_cols)

    @info "Running PAWN" n_factors=length(numeric_cols) factors=numeric_cols

    pawn_cover   = ADRIA.sensitivity.pawn(pawn_df, y_cover)
    pawn_juv     = ADRIA.sensitivity.pawn(pawn_df, y_juv)
    pawn_fish    = ADRIA.sensitivity.pawn(pawn_df, y_fish)
    pawn_tourism = ADRIA.sensitivity.pawn(pawn_df, y_tourism)
    pawn_even    = ADRIA.sensitivity.pawn(pawn_df, y_even)
    pawn_div     = ADRIA.sensitivity.pawn(pawn_df, y_div)
    pawn_shelter = ADRIA.sensitivity.pawn(pawn_df, y_shelter)
    pawn_biodiv  = ADRIA.sensitivity.pawn(pawn_df, y_biodiv)
    pawn_taxa    = ADRIA.sensitivity.pawn(pawn_df, y_taxa)

    # ── 10. PAWN figures ──────────────────────────────────────────────────────
    figures_dir = joinpath(fpath, "figures")
    mkpath(figures_dir)

    for (oname, pr) in [
        ("cover",   pawn_cover),
        ("juv",     pawn_juv),
        ("fish",    pawn_fish),
        ("tourism", pawn_tourism),
        ("even",    pawn_even),
        ("div",     pawn_div),
        ("shelter", pawn_shelter),
        ("biodiv",  pawn_biodiv),
        ("taxa",    pawn_taxa),
    ]
        fig = ADRIA.viz.pawn(pr)
        save(joinpath(figures_dir, "pawn_diff_$(oname).png"), fig; px_per_unit=3)
    end

    @info "PAWN figures saved to $figures_dir"

    return (;
        input_table,
        pawn_df,
        summary_df,
        y_cover, y_juv, y_fish, y_tourism,
        y_even, y_div, y_shelter, y_biodiv, y_taxa,
        pawn_cover, pawn_juv, pawn_fish, pawn_tourism,
        pawn_even, pawn_div, pawn_shelter, pawn_biodiv, pawn_taxa,
        figures_dir,
        out_dir,
    )
end
