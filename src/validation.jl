export ValidationResult, check_duplicates, check_ranges, check_completeness, validate_inventory

"""
    ValidationResult

Result of inventory validation checks.

# Fields
- `valid::Bool`: Whether the inventory passed all checks (no errors).
- `warnings::Vector{String}`: Warning messages (non-fatal issues).
- `errors::Vector{String}`: Error messages (fatal issues).
- `n_duplicates::Int`: Number of duplicate record groups found.
- `n_range_issues::Int`: Number of records with out-of-range values.
- `n_missing_fields::Int`: Number of records with missing required fields.
"""
struct ValidationResult
    valid::Bool
    warnings::Vector{String}
    errors::Vector{String}
    n_duplicates::Int
    n_range_issues::Int
    n_missing_fields::Int
end

"""
    check_duplicates(emissions::DataFrame) -> DataFrame

Find duplicate records in an emissions inventory based on `(FIPS, SCC, POLID)` groups.

Returns a DataFrame containing only the duplicate records, with an additional `:dup_count`
column indicating how many records share the same key. Returns an empty DataFrame if no
duplicates are found.
"""
function check_duplicates(emissions::DataFrame)
    if nrow(emissions) == 0 || !all(hasproperty.(Ref(emissions), [:FIPS, :SCC, :POLID]))
        return DataFrame(FIPS = String[], SCC = String[], POLID = String[], dup_count = Int[])
    end

    gdf = groupby(emissions, [:FIPS, :SCC, :POLID])
    counts = combine(gdf, nrow => :dup_count)
    dups = filter(r -> r.dup_count > 1, counts)
    return dups
end

"""
    check_ranges(emissions::DataFrame; max_value::Float64=1e12) -> DataFrame

Flag records with negative, zero, or extreme `ANN_VALUE` values.

Returns a DataFrame of flagged records with an additional `:range_issue` column
describing the problem (`:negative`, `:zero`, or `:extreme`).
"""
function check_ranges(emissions::DataFrame; max_value::Float64 = 1e12)
    if nrow(emissions) == 0 || !hasproperty(emissions, :ANN_VALUE)
        return DataFrame(
            FIPS = String[], SCC = String[], POLID = String[],
            ANN_VALUE = Float64[], range_issue = Symbol[]
        )
    end

    issues = DataFrame(
        FIPS = String[], SCC = String[], POLID = String[],
        ANN_VALUE = Float64[], range_issue = Symbol[]
    )

    for row in eachrow(emissions)
        val = Float64(ustrip(row.ANN_VALUE))
        issue = if val < 0
            :negative
        elseif val == 0
            :zero
        elseif abs(val) > max_value
            :extreme
        else
            nothing
        end
        if issue !== nothing
            push!(issues, (
                FIPS = string(row.FIPS),
                SCC = string(row.SCC),
                POLID = string(row.POLID),
                ANN_VALUE = val,
                range_issue = issue
            ))
        end
    end
    return issues
end

"""
    check_completeness(emissions::DataFrame;
        required::Vector{Symbol}=[:FIPS, :SCC, :POLID, :ANN_VALUE]) -> DataFrame

Flag records with missing or empty values in required fields.

Returns a DataFrame of flagged records with an additional `:missing_fields` column
listing which required fields are missing or empty.
"""
function check_completeness(emissions::DataFrame;
        required::Vector{Symbol} = [:FIPS, :SCC, :POLID, :ANN_VALUE])
    if nrow(emissions) == 0
        return DataFrame(row_index = Int[], missing_fields = Vector{Symbol}[])
    end

    issues = DataFrame(row_index = Int[], missing_fields = Vector{Symbol}[])

    for (i, row) in enumerate(eachrow(emissions))
        missing_cols = Symbol[]
        for col in required
            if !hasproperty(emissions, col)
                push!(missing_cols, col)
            elseif ismissing(row[col])
                push!(missing_cols, col)
            elseif row[col] isa AbstractString && isempty(strip(string(row[col])))
                push!(missing_cols, col)
            end
        end
        if !isempty(missing_cols)
            push!(issues, (row_index = i, missing_fields = missing_cols))
        end
    end
    return issues
end

"""
    validate_inventory(emissions::DataFrame) -> ValidationResult

Run all validation checks on an emissions inventory DataFrame.

Calls [`check_duplicates`](@ref), [`check_ranges`](@ref), and
[`check_completeness`](@ref), then aggregates results into a
[`ValidationResult`](@ref).

# Returns
A `ValidationResult` with `valid=true` if there are no errors
(missing required fields or negative values). Duplicates and zero/extreme
values are reported as warnings.
"""
function validate_inventory(emissions::DataFrame)
    warnings = String[]
    errors = String[]

    # Check duplicates
    dups = check_duplicates(emissions)
    n_dups = nrow(dups)
    if n_dups > 0
        push!(warnings, "Found $n_dups duplicate (FIPS, SCC, POLID) groups")
    end

    # Check ranges
    range_issues = check_ranges(emissions)
    n_range = nrow(range_issues)
    n_negative = count(r -> r.range_issue == :negative, eachrow(range_issues))
    n_zero = count(r -> r.range_issue == :zero, eachrow(range_issues))
    n_extreme = count(r -> r.range_issue == :extreme, eachrow(range_issues))
    if n_negative > 0
        push!(errors, "Found $n_negative records with negative ANN_VALUE")
    end
    if n_zero > 0
        push!(warnings, "Found $n_zero records with zero ANN_VALUE")
    end
    if n_extreme > 0
        push!(warnings, "Found $n_extreme records with extreme ANN_VALUE")
    end

    # Check completeness
    missing_fields = check_completeness(emissions)
    n_missing = nrow(missing_fields)
    if n_missing > 0
        push!(errors, "Found $n_missing records with missing required fields")
    end

    valid = isempty(errors)
    return ValidationResult(valid, warnings, errors, n_dups, n_range, n_missing)
end
