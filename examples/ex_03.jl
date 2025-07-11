using ControlPlots

function get_properties_along_curve(elevation_max, azimuth_max, s)
        # Elevation and azimuth as function of normalized arc length.
        theta = elevation_max * sin(4 * pi * s)  # [rad]
        phi = azimuth_max * sin(2 * pi * s)  # [rad]
        # Derivatives wrt normalized arc length.
        dtheta_ds = 4 * pi * elevation_max * cos(4 * pi * s)  # [-]
        dphi_ds = 2 * pi * azimuth_max * cos(2 * pi * s)  # [-]

        chi = atan(dphi_ds, -dtheta_ds)

        return theta, phi, chi, dtheta_ds, dphi_ds

end



elevation_max = 5 * pi / 180
azimuth_max = 30 * pi / 180

theta = 0.
phi = 0.
chi = 0.
dtheta_ds = 0.
dphi_ds = 0.

for s in LinRange(0, 2, 201)
    theta, phi, chi, dtheta_ds, dphi_ds .= get_properties_along_curve(elevation_max, azimuth_max, s)
    
    #plt.scatter(phi, theta)
    #plt.scatter(s, theta, color="C1")
    #plt.scatter(s, phi, color="C2")
    #plt.scatter(s, chi, color="C3")

end

plt.show()