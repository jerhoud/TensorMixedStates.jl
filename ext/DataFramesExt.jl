module DataFramesExt

import TensorMixedStates: data_to_frame
using DataFrames

function data_to_frame(data::Dict)
    dfs = [DataFrame("time" => identity.(val["times"]), key => identity.(val["data"])) for (key, val) in data]
    if length(dfs) == 1
      dfs[1]
    else
      outerjoin(dfs...; on = :time)
    end
end

end