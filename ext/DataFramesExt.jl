module DataFramesExt

import TensorMixedStates: data_to_frame
using TensorMixedStates: DataTable, DataRow
using DataFrames

# the rows of a Data, a slice of them or a vector collected from them, see DataTable
data_to_frame(rows::AbstractVector{DataRow}) = DataFrame(DataTable(collect(rows)))

end
