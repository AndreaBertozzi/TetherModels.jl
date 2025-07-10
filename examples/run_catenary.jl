using ControlPlots, TetherModels, StaticArrays, LinearAlgebra, KiteUtils, AtmosphericModels

# Set wind speed and direction
v_wind_gnd = 0.0
wind_dir = pi/4

# Load settings using KiteUtils
set_data_path("data")
settings = load_settings("system.yaml")

# Calculate kite position and velocity from settings
d = settings.kite_distance
ct = cos(settings.elevation)
st = sin(settings.elevation)
cp = cos(settings.azimuth)
sp = sin(settings.azimuth)

kite_pos = MVector{3}([d * ct * cp,
                    d * ct * sp,                    
                    d * st])

kite_vel = MVector{3}([0., 0., 0.])

# Change some default settings
settings.segments = 20

# Create atmospheric model using AtmosphericModels
am = AtmosphericModel()

# Define tether object and initialize with analytic catenary shape
tether = Tether(settings, am)
init_tether!(tether)

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
plt.show()

