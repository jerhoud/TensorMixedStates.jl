module DataFramesExt

import TensorMixedStates: data_to_frame
using DataFrames

# one row per call of `output`, the values measured together: the time alone repeats over the
# sweeps of dmrg, a circuit, or once it is set back, and joining on it paired values that were
# never measured together
function data_to_frame(data::Dict)
    # nothing gathered: outerjoin of no table failed
    if isempty(data)
        return DataFrame(time = [])
    end
    # a measurement named time or event keeps its column, renamed time_1 or event_1. In the
    # order of their names, which the order of the dictionary changed from one process to the
    # next
    dfs = [DataFrame("event" => identity.(val["events"]), "time" => identity.(val["times"]),
                     key => identity.(val["data"]); makeunique = true)
           for (key, val) in sort!(collect(data); by = first)]
    df = length(dfs) == 1 ? dfs[1] : outerjoin(dfs...; on = [:event, :time], makeunique = true)
    sort!(df, :event)
    return select!(df, Not(:event))
end

end