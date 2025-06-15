using ControlPlots, TetherModels, StaticArrays, LinearAlgebra, KiteUtils

# Set initial conditions
kite_pos = MVector{3}([100.0, 100, 800])
tether_length = norm(kite_pos)*1.05

set_data_path("data")
settings = load_settings("system.yaml")
settings.l_tether = norm(kite_pos)*1.05
settings.segments = 100

state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = init_quasistatic(kite_pos, tether_length, segments = 10, set=settings)

settings = tether_set_from_set(settings)
tether = Tether(settings, state_vec = state_vec)

v_wind_gnd = 0.0
wind_dir = 0.0

step!(tether, kite_pos, kite_vel, v_wind_gnd, wind_dir; prn=false)

state_vec = tether.state_vec
tether_pos = tether.tether_pos


x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])


tether_pos = hcat(tether_pos, [0; 0; 0])
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
#plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.zlim(0, 800)

plt.legend(["Tether", "Origin", "Kite"])
plt.show()

