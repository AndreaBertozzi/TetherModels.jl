using ControlPlots, TetherModels, StaticArrays, LinearAlgebra, KiteUtils

# Set initial conditions
kite_pos = MVector{3}([100.0, 100, 800])
tether_length = norm(kite_pos)*1.05

set_data_path("data")
settings = load_settings("system.yaml")

state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = init_quasistatic(kite_pos, tether_length, segments = 22, settings=settings)
state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)


x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])


tether_pos = hcat(p0, tether_pos, [0; 0; 0])
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.zlim(0, 800)

plt.legend(["Tether", "Origin", "Kite"])
plt.show()

