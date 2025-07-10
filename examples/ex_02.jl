using ControlPlots, TetherModels, StaticArrays, LinearAlgebra, KiteUtils, AtmosphericModels
"""
Example 1: static calculation of a slack tether. The tether is first initialized
solving the analytic catenary equation. The quasi-static model is then runned. Both solutions
are plotted superimposed one on the other. The expected resuls is a plot with the blue line
(analytic catenary solution) and the orange line (quasi-static numerical solution) very
close to each other. 
"""
# Define helper function
    """Calculate kite position and velocity from settings"""
function kite_pos_from_settings(set::Settings)    
    d = set.kite_distance
    ct = cos(set.elevation)
    st = sin(set.elevation)
    cp = cos(set.azimuth)
    sp = sin(set.azimuth)

    kite_pos = MVector{3}([d * ct * cp,
                           d * ct * sp,                    
                           d * st])
end

# Load system settings
set_data_path("data")
set = load_settings("system.yaml")

# Change kite distance to have slack tether 
set.segments = 15
set.azimuth = 0.3
set.kite_distance = 49.8

# Get kite position from settings
kite_pos = kite_pos_from_settings(set)
kite_vel = MVector{3}([0., 0., 0.])

# Create atmospheric model
am = AtmosphericModel(set)
# Set zero wind
v_wind_gnd= 0.
wind_dir = 0.

# Define tether object
tether = Tether(set, am)

# Initalize tether object. In this example, kite_distance < l_tether, 
# so the tether is initialized under straight tether approximation
init_tether!(tether)



# Plot initialised straight tether trajectory and nodes positions
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.plot(tether.tether_pos[1,:], tether.tether_pos[2,:], tether.tether_pos[3,:], marker = "+", linewidth = 2)

# Run model iteration 
step!(tether, kite_pos, kite_vel, v_wind_gnd, wind_dir; prn=false)

# Plot model after iteration
plt.plot(tether.tether_pos[1,:], tether.tether_pos[2,:], tether.tether_pos[3,:], marker = "x", linestyle = "--")
plt.scatter3D(0, 0, 0,  s = 50, marker = "s", c = "C7")
plt.scatter3D(kite_pos[1], kite_pos[2], kite_pos[3],  s = 50, marker = "^", c = "C2")

plt.legend(["Analytic catenary initialization", "Resulting model iteration", "Ground station", "Kite"])
#plt.xlim(0, 50)
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.show()

