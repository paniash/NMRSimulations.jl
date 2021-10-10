using DifferentialEquations
using Plots

##########################################
#%%
function magnetization!(du, u, p, t)
    du[1] = -u[1]/p[2] + (p[3]) * u[2]
    du[2] = -p[3] * u[1] - u[2]/p[2] - p[4]*u[3]
    du[3] = p[4]*u[2] - u[3]/p[1] + 1/p[1]
end

#%%
T1 = 1; T2 = 0.3
ω1 = 2*π*10
# ω1 = 0
Δ = 0.5 * 2 * π
p = [T1, T2, Δ, ω1]
u0 = [0, 1, 0]
tmax = 0.5
tspan = (0.0, tmax)

#%%
# Define function to plot components of magnetization
function components_magnetization(func,u0,tspan,p)
    problem = ODEProblem(func,u0,tspan,p)
    sol = solve(problem, RK4(), saveat=0.001)
    t = sol.t
    u = sol[1,:]
    v = sol[2,:]
    w = sol[3,:]
    return t, u, v, w
end

#%% Solve for variations in T1, T2
tspan = (0.0, 5.0)
t1, u1, v1, w1 = components_magnetization(magnetization!, u0, tspan, [1, 0.3, 0.5*2*π, 2*π*10])
t2, u2, v2, w2 = components_magnetization(magnetization!, u0, tspan, [10, 10, 0.5*2*π, 2*π*10])
t3, u3, v3, w3 = components_magnetization(magnetization!, u0, tspan, [Inf, Inf, 0.5*2*π, 2*π*10])
t4, u4, v4, w4 = components_magnetization(magnetization!, [0,0,1], (0.0, 5.0), [1, 0.3, 0.5*2*π, 0])

#%% Plot u(t) vs t
plot(t1, u1, xlims=(0,2), label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(t2, u2, xlims=(0,2), label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot!(t3, u3, xlims=(0,2), label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
xlabel!("t (sec)")
ylabel!("u(t)")

plot(t4, u4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("t (s)")
ylabel!("u(t)")

#%% Plot v(t) vs t
plot(t1, v1, xlims=(0,3), label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(t2, v2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot(t3, v3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(t4, v4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("t (s)")
ylabel!("v(t)")

#%% Plot w(t) vs t
plot(t1, w1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(t2, w2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot(t3, w3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(t4, w4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("t (s)")
ylabel!("w(t)")

#%% Plot u(t) vs v(t)
plot(u1, v1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(u2, v2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot!(u3, v3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(u4, v4, label="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("u(t)")
ylabel!("v(t)")

#%% Plot u(t) vs w(t)
plot(u1, w1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(u2, w2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot!(u3, w3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(u4, w4, label="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("u(t)")
ylabel!("w(t)")
legend()

#%% Plot v(t) vs w(t)
plot(v1, w1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot!(v2, w2, label="T1 = 10s, T2 = 10 s", lw=1.5, grid=false)
plot!(v3, w3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(v4, w4, label="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("v(t)")
ylabel!("w(t)")

#%% 3D plot for u(t), v(t), w(t)
plot(u1, w1, v1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot(u2, w2, v2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot(u3, w3, v3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(u4, v4, w4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("u(t)")
ylabel!("v(t)")

#%% 3D plot for u(t), v(t), t
plot(u1, v1, t1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot(u2, v2, t2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot(u3, v3, t3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(u4, v4, t4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("u(t)")
ylabel!("v(t)")

#%% 3D plot for v(t), w(t), t
plot(v1, w1, t1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot(v2, w2, t2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot(v3, w3, t3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(t4, v4, w4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("v(t)")
ylabel!("w(t)")

#%% 3D plot for u(t), w(t), t
plot(u1, w1, t1, label="T1 = 1 s, T2 = 0.3 s", lw=1.5, grid=false)
plot(u2, w2, t2, label="T1 = 10 s, T2 = 10 s", lw=1.5, grid=false)
plot!(u3, w3, t3, label="T1 = ∞, T2 = ∞", lw=1.5, grid=false)
plot(t4, u4, v4, title="Free precession for T1 = 1s, T2 = 0.3s", lw=1.5, grid=false)
xlabel!("u(t)")
ylabel!("w(t)")

#%%
# Δ = range(0.0, 20.0, step=0.01)

# function lin_solve(delta, A, b, T2)
#     A = [1 -T2*delta; -T2*delta -1]
#     b = [0; 0]
#     return A\b[1], A\b[2]
# end

# for i=1:48
#     u = []; v = []
#     push!(u, lin_solve(Δ[i], [1 -T2*Δ[i]; -T2*Δ[i] -1], [0; 0], T2)[1])
#     push!(v, lin_solve(Δ[i], [1 -T2*Δ[i]; -T2*Δ[i] -1], [0; 0], T2)[2])
# end
