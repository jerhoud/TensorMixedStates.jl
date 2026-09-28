module DataFramesExt

import TensorMixedStates: data_to_frame
using DataFrames

# one row per call of `output`, the values measured together: the time alone repeats over the
# sweeps of dmrg, a circuit, or once it is set back, and joining on it paired values that were
# never measured together
function data_to_frame(data::Dict)
    dfs = [DataFrame("event" => identity.(val["events"]), "time" => identity.(val["times"]),
                     key => identity.(val["data"])) for (key, val) in data]
    df = length(dfs) == 1 ? dfs[1] : outerjoin(dfs...; on = [:event, :time])
    sort!(df, :event)
    return select!(df, Not(:event))
end

end