module DataFramesExt

import TensorMixedStates: DataToFrame
using DataFrames

function DataToFrame(data::Dict)
    dfs = [DataFrame("time" => identity.(val["times"]), key => identity.(val["data"])) for (key, val) in data]
    if length(dfs) == 1
      dfs[1]
    else
      outerjoin(dfs...; on = :time)
    end
end

end