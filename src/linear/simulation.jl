"""
    Linear_syst_int_obs_ode!(dz, z, p, t)

ODE function for linear interval observer dynamics.

Implements the coupled ODE system for a linear observer:
  ẋx = A*x
  ẋxᴶ = A*xᴶ + K*(y - C*xᴵ)
  ẋxᴵ = A*xᴵ + K*(y - C*xᴶ)

where:
- x is the true state (when available)
- xᴶ, xᴵ are the upper and lower estimates
- K is the observer gain
- y = C*x is the measurement

State vector layout: z = [x; xᴶ; xᴵ] (each of length n)

# Arguments
- `dz`: Derivative of state (output)
- `z`: Current state [x; xᴶ; xᴵ]
- `p`: Parameters tuple: (A, C, K, n)
  - `A::Matrix`: System matrix
  - `C::Vector`: Measurement vector
  - `K::Vector`: Observer gain
  - `n::Int`: System dimension
- `t`: Current time (unused in autonomous system)

# Note
This function is in-place: it modifies dz. Used as the ODE function for DifferentialEquations.jl.
"""
function Linear_syst_int_obs_ode!(dz, z, p, t)
    A, C, K, n = p

    @inbounds begin

        x =view(z, 1:n)
        xu = view(z, n+1:2n)
        xl = view(z, 2n+1:3n)

        dx = view(dz, 1:n)
        dxu = view(dz, n+1:2n)
        dxl = view(dz, 2n+1:3n)
    end   

    y = dot(C, x)
    ### this part allocates more or less
    dx .= A*x
    dxl .= A*xl + K*(y-dot(C, xu))
    dxu .= A*xu + K*(y-dot(C, xl))
end



