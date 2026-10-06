module DataFramesExt

import TensorMixedStates: data_to_frame
using DataFrames

# one row per call of `output`: the time alone repeats, over the sweeps of dmrg for instance,
# and joining on it would pair values never measured together
function data_to_frame(data::Dict)
    # nothing gathered: no table to join
    if isempty(data)
        return DataFrame(time = [])
    end
    # a measurement named time or event keeps its column, renamed time_1 or event_1; sorted by
    # name, the order of a dictionary depending on the order its keys were inserted in
    dfs = [DataFrame("event" => identity.(val["events"]), "time" => identity.(val["times"]),
                     key => identity.(val["data"]); makeunique = true)
           for (key, val) in sort!(collect(data); by = first)]
    df = length(dfs) == 1 ? dfs[1] : outerjoin(dfs...; on = [:event, :time], makeunique = true)
    sort!(df, :event)
    return select!(df, Not(:event))
end

end