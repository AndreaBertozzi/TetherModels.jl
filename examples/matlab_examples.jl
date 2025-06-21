"""
    get_analytic_catenary(filename)

Loads the analytic catenary curve for the 2D catenary example

# Arguments
- filename: the filename of the mat file to read

# Returns
- x_cat: x coordinates of the catenary curve
- y_cat: x coordinates of the catenary curve
"""
function get_analytic_catenary(filename)
    vars        = matread(filename)
    vars        = get(vars, "analytic_catenary", 0)
    x_cat       = vec(get(vars, "x", 0))
    y_cat       = vec(get(vars, "y", 0))
    return x_cat, y_cat
end

function transformFromOtoW(windDirection_rad,vec_O)
    M_WO = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0;  0 0 -1]
    vec_W = M_WO*vec_O
    return vec_W
end

function transformFromWtoO(windDirection_rad,vec_W)
    M_OW = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0; 0 0 -1]
    vec_O = M_OW*vec_W
    return vec_O
end


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
