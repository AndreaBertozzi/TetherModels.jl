module TetherModels
# TODO: Implement clear interface 

using LinearAlgebra, StaticArrays, ADTypes, NonlinearSolve, KiteUtils, AtmosphericModels

export Tether, step!, init_tether!

const MVec3 = MVector{3, Float64}
const SVec3 = SVector{3, Float64}
const GRAV_ACC = SVec3(0.0, 0.0, -9.81)

# Iterations: 36
# BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
#  Range (min … max):  84.169 μs …  2.388 ms  ┊ GC (min … max): 0.00% … 93.18%
#  Time  (median):     89.110 μs              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   92.513 μs ± 47.339 μs  ┊ GC (mean ± σ):  1.10% ±  2.09%

#     ▂█▇▆▂▁                                                     
#   ▁▃██████▆▆▆▅▅▅▄▅▄▄▄▃▃▄▃▃▃▃▃▃▄▃▃▃▃▃▃▃▃▃▃▃▃▂▂▂▃▂▂▂▂▂▂▂▂▂▂▂▁▁▁ ▃
#   84.2 μs         Histogram: frequency by time         108 μs <

#  Memory estimate: 41.20 KiB, allocs estimate: 971.

mutable struct Tether{N}
    set::Settings
    am::AtmosphericModel
    state_vec::MVec3
    tether_pos::MMatrix{3, N}
    wind_matrix::MMatrix{3, N}
    
    function Tether{N}(set::Settings, am::AtmosphericModel;
                       state_vec::MVec3 = MVector{3}(0.0, 0.0, 0.0)) where {N}
        tether_pos = MMatrix{3, N}(zeros(3, N))
        wind_matrix = MMatrix{3, N}(zeros(3, N))
        new(set, am, state_vec, tether_pos, wind_matrix)
    end
end

function Tether(set::Settings, am::AtmosphericModel;
                state_vec::MVec3 = MVector{3}(0.0, 0.0, 0.0))
    Tether{set.segments}(set, am; state_vec)
end

"""
    update_wind_matrix!(tether::Tether)

Updates the wind velocity components in the `tether.wind_matrix` for each segment in the tether.

The function computes the wind velocity (`v_wind`) based on the height of each tether segment and 
adjusts the horizontal wind components (`u` and `v`) based on a given wind direction (`wind_dir`). 

The wind velocity is calculated differently depending on the height of each tether segment:
- If the height is positive, `v_wind` is computed using the `calc_wind_factor` function,
    which factors in the tether's properties and wind profile.
- If the height is zero or negative, the wind velocity is assumed to be constant (`v_wind_gnd`).

# Arguments
- `tether::Tether`: A `Tether` object

# Returns
- The function updates the `tether.wind_matrix` in place.
"""

function update_wind_matrix!(tether::Tether, v_wind_gnd, wind_dir)
    am = tether.am
    profile_law = am.set.profile_law
    cos_wind_dir = cos(wind_dir)
    sin_wind_dir = sin(wind_dir)
    
    for i in 1:tether.set.segments
        # Calculate wind velocity for this segment
        height = tether.tether_pos[3, i]
        if height > 0
            v_wind = v_wind_gnd * calc_wind_factor(am, height, profile_law)
        else
            v_wind = v_wind_gnd
        end

        u = v_wind * cos_wind_dir
        v = v_wind * sin_wind_dir

        tether.wind_matrix[1, i] = u
        tether.wind_matrix[2, i] = v
    end       
end

"""
    simulate_tether(state_vec, kite_pos, kite_vel, wind_matrix, l_tether, set)

Function to determine the tether shape and forces, based on a quasi-static model.

# Arguments
- state_vec::MVector{3, Float64}: state vector (theta [rad], phi [rad], Tn [N]);  
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64}: kite position vector in wind reference frame
- kite_vel::MVector{3, Float64}: kite velocity vector in wind reference frame
- wind_matrix:: (3, segments) MMatrix{Float64} wind velocity vector in wind reference frame for each segment of the tether
- l_tether:: Float64: tether length
- set:: TetherSettings or Settings: struct containing environmental and tether parameters: see [TetherSettings](@ref)
                or Settings in KiteUtils.jl

# Returns
- state_vec::MVector{3, Float64}: state vector (theta [rad], phi [rad], Tn [N]);  
  tether orientation and tension at ground station 
- tether_pos::Matrix{Float64}: x,y,z - coordinates of the tether nodes
- force_gnd::Float64: Line tension at the ground station
- force_kite::Vector{Float64}: force from the kite to the end of tether
- p0::Vector{Float64}:  x,y,z - coordinates of the kite-tether attachment
"""
function step!(tether, kite_pos, kite_vel, v_wind_gnd, wind_dir; prn=false)
    
    segments = tether.set.segments
    buffers= [MMatrix{3, segments}(zeros(3, segments)), MMatrix{3, segments}(zeros(3, segments)),
              MMatrix{3, segments}(zeros(3, segments)), MMatrix{3, segments}(zeros(3, segments)),
              MMatrix{3, segments}(zeros(3, segments))]
    
    update_wind_matrix!(tether, v_wind_gnd, wind_dir)
    # Pack parameters in param named tuple - false sets res! for in-place solution
    param = (tether = tether, kite_pos=kite_pos, kite_vel=kite_vel,
                 buffers=buffers, return_result=false)
    # Define the nonlinear problem
    prob = NonlinearProblem(res!, tether.state_vec, param)
    # Solve the problem with TrustRegion method
    sol = solve(prob, TrustRegion(autodiff=AutoFiniteDiff()), show_trace=Val(false)) 
    
    if prn
        iterations = sol.stats.nsteps  # Field name may vary; verify with `propertynames(sol)`
        println("Iterations: ", iterations)
    end
    tether.state_vec = sol.u
    # Set the return_result to true so that res! returns outputs
    param = (; param..., return_result=true)
    res = MVector(0., 0., 0.)
    res, force_kite, tether.tether_pos, p0 = res!(res, tether.state_vec, param)
end


"""
    res!(res, state_vec, param)

Calculates difference between tether end and kite given tether ground segment orientation 
and magnitude.

# Arguments
- res::Vector{Float64} difference between tether end and kite segment
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N]);
  tether orientation and tension at ground station
- par:: 7-elements tuple:
    - kite_pos::MVector{3, Float64} kite position vector in wind reference frame
    - kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
    - wind_matrix::MMatrix{Float64} wind velocity vector in wind reference frame for each segment of the tether
    - l_tether: tether length
    - set:: TetherSettings struct containing environmental and tether parameters: see [TetherSettings](@ref)
    - buffers:: (5, ) Vector{Matrix{Float64}}  Vector of (3, segments) Matrix{Float64} empty matrices for preallocation
    - segments:: number of tether segments
    - return_result:: Boolean to determine use for in-place optimization or for calculating returns

# Returns (if return_result==true)
- res::Vector{Float64} difference between tether end and kite segment
- T0::Vector{Float64} force from the kite to the end of tether
- pj:: (3, segments) Matrix{Float64} x,y,z - coordinates of the tether nodes
- p0::Vector{Float64}  x,y,z - coordinates of the kite-tether attachment

# Example usage
state_vec = rand(3,)
kite_pos = [100, 100, 300] 
kite_vel = [0, 0, 0]
wind_matrix = rand(3,15)
l_tether = 500
set = TetherSettings(1.225, [0, 0, -9.806], 0.9, 4, 0.85, 500000)
res!(res, state_vec, kite_pos, kite_vel, wind_matrix, l_tether, set)
"""
function res!(res, state_vec, param)
    tether, kite_pos, kite_vel, buffers, return_result = param
    set = tether.set
    segments = tether.set.segments

    g = abs(GRAV_ACC[3])

    Ls = set.l_tether / (segments + 1)
    A = π/4 * (set.d_tether/1000)^2
    mj = set.rho_tether * Ls * A
    E = set.c_spring / A

    # Preallocate arrays
    FT = buffers[1]
    Fd = buffers[2]
    pj = buffers[3]
    vj = buffers[4]
    aj = buffers[5]

    # Unpack state variables (elevation, azimuth)
    θ, φ, Tn = state_vec[1], state_vec[2], state_vec[3]

    # Precompute common values
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)
    norm_p = norm(kite_pos)
    p_unit = kite_pos ./ norm_p
    v_parallel = dot(kite_vel, p_unit)
    
    # First element calculations
    FT[1, segments] = Tn * cosθ * cosφ # cos(elevation)cos(azimuth)
    FT[2, segments] = Tn * cosθ * sinφ # cos(elevation)sin(azimuth)
    FT[3, segments] = Tn * sinθ        # sin(elevation)

    pj[1, segments] = Ls * cosθ * cosφ
    pj[2, segments] = Ls * cosθ * sinφ
    pj[3, segments] = Ls * sinθ

    # Velocity and acceleration calculations
    ω = cross(kite_pos / norm_p^2, kite_vel)
    a = cross(ω, SVec3(pj[:, segments]))         
    b = cross(ω, cross(ω, SVec3(pj[:, segments])))
    vj[:, segments] .= v_parallel * p_unit + a
    aj[:, segments] .= b

    # Drag calculation for first element
    v_a_p1 = vj[1, segments] - tether.wind_matrix[1, segments]
    v_a_p2 = vj[2, segments] - tether.wind_matrix[2, segments]
    v_a_p3 = vj[3, segments] - tether.wind_matrix[3, segments]

    if all(x -> abs(x) < 1e-3, (v_a_p1, v_a_p2, v_a_p3))
        Fd[:, segments] .= 0.0
    else
        dir1, dir2, dir3 = pj[1, segments]/Ls, pj[2, segments]/Ls,
                     pj[3, segments]/Ls

        v_dot_dir = v_a_p1*dir1 + v_a_p2*dir2 + v_a_p3*dir3
        v_a_p_t1 = v_dot_dir * dir1
        v_a_p_t2 = v_dot_dir * dir2
        v_a_p_t3 = v_dot_dir * dir3

        v_a_p_n1 = v_a_p1 - v_a_p_t1
        v_a_p_n2 = v_a_p2 - v_a_p_t2
        v_a_p_n3 = v_a_p3 - v_a_p_t3

        norm_v_a_p_n = sqrt(v_a_p_n1^2 + v_a_p_n2^2 + v_a_p_n3^2)
        drag_coeff = -0.5 * set.rho_0 * Ls * set.d_tether * set.cd_tether
        coeff = drag_coeff * norm_v_a_p_n

        Fd[1, segments] = coeff * v_a_p_n1
        Fd[2, segments] = coeff * v_a_p_n2
        Fd[3, segments] = coeff * v_a_p_n3
    end

    # Process other segments
    @inbounds for ii in segments:-1:2
        # Tension force calculations
        if ii == segments
            mj_total = 1.5mj
            g_term = mj_total * g
        else
            mj_total = mj
            g_term = mj * g
        end

        FT[:, ii-1] .= mj_total * aj[:, ii] + FT[:, ii] - Fd[:, ii]
        FT[3, ii-1] += g_term

        # Position calculations
        ft_norm = sqrt(FT[1, ii-1]^2 + FT[2, ii-1]^2 + FT[3, ii-1]^2)
        l_i_1 = (ft_norm/(E*A) + 1) * Ls
        ft_dir = FT[1, ii-1]/ft_norm, FT[2, ii-1]/ft_norm, FT[3, ii-1]/ft_norm

        pj[1, ii-1] = pj[1, ii] + l_i_1 * ft_dir[1]
        pj[2, ii-1] = pj[2, ii] + l_i_1 * ft_dir[2]
        pj[3, ii-1] = pj[3, ii] + l_i_1 * ft_dir[3]

        # Velocity and acceleration
        a = cross(ω, SVec3(pj[:, ii-1]))           
        b = cross(ω, cross(ω, SVec3(pj[:, ii-1])))
        vj[:, ii-1] .= v_parallel * p_unit + a
        aj[:, ii-1] .= b

        # Drag calculations
        v_a_p1 = vj[1, ii] - tether.wind_matrix[1, ii]
        v_a_p2 = vj[2, ii] - tether.wind_matrix[2, ii]
        v_a_p3 = vj[3, ii] - tether.wind_matrix[3, ii]

        if all(x -> abs(x) < 1e-3, (v_a_p1, v_a_p2, v_a_p3))
            Fd[:, ii-1] .= 0.0
        else
            dx = pj[1, ii-1] - pj[1, ii]
            dy = pj[2, ii-1] - pj[2, ii]
            dz = pj[3, ii-1] - pj[3, ii]
            segment_norm = sqrt(dx^2 + dy^2 + dz^2)
            dir1 = dx/segment_norm
            dir2 = dy/segment_norm
            dir3 = dz/segment_norm

            v_dot_dir = v_a_p1*dir1 + v_a_p2*dir2 + v_a_p3*dir3
            v_a_p_t1 = v_dot_dir * dir1
            v_a_p_t2 = v_dot_dir * dir2
            v_a_p_t3 = v_dot_dir * dir3

            v_a_p_n1 = v_a_p1 - v_a_p_t1
            v_a_p_n2 = v_a_p2 - v_a_p_t2
            v_a_p_n3 = v_a_p3 - v_a_p_t3

            norm_v_a_p_n = sqrt(v_a_p_n1^2 + v_a_p_n2^2 + v_a_p_n3^2)
            
            height = pj[3, ii-1]
            if height > 0
                rho_at_height = calc_rho(tether.am, height) 
            else
                rho_at_height = set.rho_0
            end

            drag_coeff = -0.5 * rho_at_height * Ls * set.d_tether * set.cd_tether        
            coeff = drag_coeff * norm_v_a_p_n

            Fd[1, ii-1] = coeff * v_a_p_n1
            Fd[2, ii-1] = coeff * v_a_p_n2
            Fd[3, ii-1] = coeff * v_a_p_n3
        end
    end

    # Final ground connection calculations
    T0_1 = 1.5mj*aj[1,1] + FT[1,1] - Fd[1,1]
    T0_2 = 1.5mj*aj[2,1] + FT[2,1] - Fd[2,1]
    T0_3 = 1.5mj*aj[3,1] + FT[3,1] - Fd[3,1] + 1.5mj*g
    T0_norm = sqrt(T0_1^2 + T0_2^2 + T0_3^2)
    
    l_i_1 = (T0_norm/(E*A) + 1) * Ls
    T0_dir1 = T0_1/T0_norm
    T0_dir2 = T0_2/T0_norm
    T0_dir3 = T0_3/T0_norm

    p0 = MVector(pj[1,1] + l_i_1*T0_dir1, 
                 pj[2,1] + l_i_1*T0_dir2,
                 pj[3,1] + l_i_1*T0_dir3)

    res .= kite_pos - p0
    if return_result
        return res, MVector(T0_1, T0_2, T0_3), pj, p0
    else
        nothing
    end
end

"""
    init_quasistatic(kite_pos, l_tether; kite_vel = nothing, segments = nothing, wind_matrix = nothing, set = nothing)

Initialize the quasi-static tether model providing an initial guess for the state vector based on the numerical solution of the catenary equation

# Arguments
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- l_tether: Float64 tether length
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- segments::Int number of tether segments
- wind_matrix::MMatrix{3, segments, Float64} wind velocity vector in wind reference frame for each segment of the tether
- set::TetherSettings or Settings: struct containing environmental and tether parameters: see [TetherSettings](@ref)
                                        or Settings in KiteUtils.jl

# Returns
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N])  
  tether orientation and tension at ground station
"""
function init_tether!(tether, kite_pos; kite_vel=nothing)
    # Assert the correct types
    @assert isa(tether, Tether) "tether should be a Tether object!"
    @assert isa(kite_pos, MVector{3}) "kite_pos must be a MVector of size (3,1)"

    # Set kite velocity if not provided
    if isnothing(kite_vel)
        kite_vel = MVector{3}([0.0, 0.0, 0.0])
    else
        @assert isa(kite_vel, MVector{3}) "kite_vel must be a MVector of size (3,1)"
    end

    # Extract tether length and kite position components
    tether_length = tether.set.l_tether
    horizontal_pos = kite_pos[1:2]  # horizontal vector (x, y)
    horizontal_distance = norm(horizontal_pos)  # horizontal distance from origin
    vertical_pos = kite_pos[3]       # vertical position (z)

    # Define the nonlinear solver function for catenary
    function catenary_residual!(residual, coeff, params)
        tether_length, vertical_pos, horizontal_distance = params
        residual[] = sqrt(tether_length^2 - vertical_pos^2) - 
                        (2 * sinh(coeff[] * horizontal_distance / 2) / coeff[])
    end

    initial_guess = [0.1]  # Initial guess for the catenary coefficient as a scalar in an array

    # Solve the nonlinear problem for the catenary coefficient using Newton's method
    params = (tether_length, vertical_pos, horizontal_distance)
    prob = NonlinearProblem(catenary_residual!, initial_guess, params)
    catenary_coefficient = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()), show_trace=Val(false)) 
    catenary_coeff_val = catenary_coefficient[]  # Extract scalar value from the solution

    # Define the catenary solution based on the horizontal position and coefficient
    horizontal_positions = LinRange(0, horizontal_distance, tether.set.segments)  
    azimuth_angle = atan(horizontal_pos[2], horizontal_pos[1]) 

    # Calculate horizontal projection of the catenary curve (XY positions)
    XY_positions = [sin(azimuth_angle) * horizontal_positions'; cos(azimuth_angle) * horizontal_positions']

    # Calculate x_min and bias for the catenary curve's vertical displacement
    x_min = -(1 / 2) * (log((tether_length + vertical_pos) / (tether_length - vertical_pos)) /
                             catenary_coeff_val - horizontal_distance)
    vertical_bias = -cosh(-x_min * catenary_coeff_val) / catenary_coeff_val    

    # Assign the x, y, and z positions to the tether positions
    tether.tether_pos[1, :] .= XY_positions[1, :]  # x positions
    tether.tether_pos[2, :] .= XY_positions[2, :]  # y positions
    tether.tether_pos[3, :] .= cosh.((horizontal_positions .- x_min) .* catenary_coeff_val) ./
                                         catenary_coeff_val .+ vertical_bias  # z positions

    # Calculate the azimuth angle (phi) based on kite position
    azimuth_phi_init = atan(kite_pos[2], kite_pos[1])

    # Calculate the elevation angle (theta) from the tether's position (z, x, y)
    elevation_theta_init = atan(tether.tether_pos[3, 2],
                                sqrt(tether.tether_pos[1, 2]^2 + tether.tether_pos[2, 2]^2))

    # Calculate initial tension using the spring constant
    initial_tension = 2e-4 * tether.set.c_spring

    # Update tether state vector with initial angles and tension
    tether.state_vec = MVector{3}([elevation_theta_init, azimuth_phi_init, initial_tension]) 
end