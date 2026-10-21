module DataFramesExt

import TensorMixedStates: data_to_frame
using TensorMixedStates: DataTable, DataRow
using DataFrames

data_to_frame(rows::DataTable) = DataFrame(rows)

# rows of pairs, as a Data holds them, their columns named as values of no measurement of
# their own
data_to_frame(rows::AbstractVector) = DataFrame(DataTable(DataRow[ DataRow(r) for r in rows ]))

end
