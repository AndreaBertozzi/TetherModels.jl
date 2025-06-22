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

# Define tether object and initialize with analytic catenary shape
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


