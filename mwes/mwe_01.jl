using TetherModels
using KiteUtils

set_data_path("data")
set = load_settings("system.yaml")
println(set.d_tether)

ts = TetherModels.TetherSettings()
ts.cd_tether = set.cd_tether
ts.d_tether = set.d_tether
ts.rho_tether = set.rho_tether


# TODO calculate rho as function of height
# TODO calculate c_spring as function of d_tether, rho_tether and segment length
