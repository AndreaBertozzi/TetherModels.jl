using TetherModels
using KiteUtils

set_data_path("data")
set = load_settings("system.yaml")
println(set.d_tether)
