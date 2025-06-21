using ControlPlots, TetherModels, StaticArrays, LinearAlgebra, KiteUtils, AtmosphericModels

# Set example kite_pos and kite_vel
kite_pos = MVector{3}([100.0, 100.0, 300.0])
kite_vel = MVector{3}([0.5, 0., 0])

# Set wind speed and direction
v_wind_gnd = 1.0
wind_dir = pi/4

# Load settings using KiteUtils
set_data_path("data")
settings = load_settings("system.yaml")

# Change some default settings
settings.l_tether = norm(kite_pos)*1.05
settings.segments = 20

# Create atmospheric model using AtmosphericModels
am = AtmosphericModel()

# Define tether object and initialise with analytic catenary shape
tether = Tether(settings, am)
init_tether!(tether, kite_pos; kite_vel)

# Plot initial condition of catenary
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.plot3D(tether.tether_pos[1,:], tether.tether_pos[2,:], tether.tether_pos[3,:], marker = "+")

# Run model iteration 
step!(tether, kite_pos, kite_vel, v_wind_gnd, wind_dir; prn=false)

# Plot model after iteration
plt.plot3D(tether.tether_pos[1,:], tether.tether_pos[2,:], tether.tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.zlim(0, 350)
plt.show()


"""
    get_initial_conditions(filename)

Loads the initialization data for the basic examples and tests

# Arguments
- filename: the filename of the mat file to read

# Returns
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N])  
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- wind_matrix::MMatrix{3, segments, Float64} wind velocity vector in wind reference frame for each segment of the tether
- l_tether: Float64 tether length
- set::TetherSettings struct containing environmental and tether parameters: see [TetherSettings](@ref)
"""
function get_initial_conditions(filename)
    vars = matread(filename) 
    state_vec = MVector{3}(vec(get(vars,"stateVec", 0)))
    kite_pos = MVector{3}(vec(get(vars,"kitePos", 0)))
    kite_vel = MVector{3}(vec(get(vars,"kiteVel", 0)))
    wind_matrix = get(vars,"windVel", 0)
    l_tether = get(vars,"tetherLength", 0)

    ENVMT = get(vars,"ENVMT", 0) 
    rho_air = get(ENVMT, "rhos", 0) 
    g_earth = [0; 0; -abs(get(ENVMT, "g", 0))]      # in this way g_earth is a vector [0; 0; -9.81]

    T = get(vars,"T", 0)
    cd_tether = get(T, "CD_tether", 0) 
    d_tether = get(T, "d_tether", 0)*1000           # tether diameter                  [mm]
    rho_tether = get(T, "rho_t", 0) 
    E = get(T, "E", 0) 
    A = get(T, "A", 0)
    c_spring = E*A 

    set = TetherSettings(rho_air, cd_tether, d_tether, rho_tether, c_spring)

    return state_vec, kite_pos, kite_vel, wind_matrix, l_tether, set
end
