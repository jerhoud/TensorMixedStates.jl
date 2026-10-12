module DataFramesExt

import TensorMixedStates: data_to_frame
using DataFrames

# one row per call of output, the measurements in the order they first come, a list written
# later to the same Data adding its columns, missing in the rows before
function data_to_frame(rows::AbstractVector)
    # nothing gathered: no column but the times
    if isempty(rows)
        return DataFrame(time = [])
    end
    names = unique(first(p) for r in rows for p in r)
    found = [ Dict(r) for r in rows ]
    return DataFrame([ n => map(identity, [ get(f, n, missing) for f in found ]) for n in names ])
end

end
