using DifferentialEquations
using FFTW
using Plots

# Define the coupled ODEs
function bloch!(du, u, p, t)
    du[1] = 0.5 * (p[2]*u[2] + p[3] * cos(p[1]*t) * u[4])
    du[2] = -0.5 * (p[2] * u[1] + p[3] * cos(p[1]*t) * u[3])
    du[3] = -0.5 * (p[2] * u[4] - p[3] * cos(p[1]*t) * u[2])
    du[4] = 0.5 * (p[2] * u[3] - p[3] * cos(p[1]*t) * u[1])
end

function prob_orientation(func, u0, tspan, p)
    problem = ODEProblem(func,u0,tspan,p)
    sol = solve(problem,RK4(),saveat=0.1)
    probd = sol[3,:].^2 + sol[4,:].^2
    probu = sol[1,:].^2 + sol[2,:].^2
    return sol.t, probd, probu
end

#%%
u0 = [1.0;0;0;0]
tspan = (0.0,200.0)
p = [1.0;0.95;0.1]

plot(prob_orientation(bloch!, u0, tspan, p)[1], prob_orientation(bloch!, u0, tspan, p)[2], label="Spin down", lw=1.5, grid=false)
plot!(prob_orientation(bloch!, u0, tspan, p)[1], prob_orientation(bloch!, u0, tspan, p)[3], label="Spin up", lw=1.5, grid=false)
xlabel!("t (ns)")
ylabel!("Probability of spin orientation")

#%%
## Plotting flipping probability
N = 1000
ω0 = 1.0; ω1 = 0.1
tmax = 2.1 * pi / ω1
tspan = (0.0,tmax)
omega_array = range(0.9, 1.1, length=N)
prob_array1 = zeros(N)

for i = 1:N
    ω = omega_array[i]
    p = [ω, ω0, ω1]
    p_b = prob_orientation(bloch!, u0, tspan, p)[2]
    prob_array1[i] = maximum(p_b)
end

# Repeating the above but for a different value of ω1
#%%
N = 1000
ω0 = 1.0; ω1 = 0.01
tmax = 2.1 * pi / ω1
tspan = (0.0,tmax)
omega_array = range(0.9, 1.1, length=N)
prob_array2 = zeros(N)

for i = 1:N
    ω = omega_array[i]
    param = [ω, ω0, ω1]
    p_b = prob_orientation(bloch!, u0, tspan, param)[2]
    prob_array2[i] = maximum(p_b)
end

# Repeating the above but for a different value of ω1
#%%
N = 1000
ω0 = 1.0; ω1 = 0.05
tmax = 2.1 * pi / ω1
tspan = (0.0,tmax)
omega_array = range(0.9, 1.1, length=N)
prob_array3 = zeros(N)

for i = 1:N
    ω = omega_array[i]
    param = [ω, ω0, ω1]
    p_b = prob_orientation(bloch!, u0, tspan, param)[2]
    prob_array3[i] = maximum(p_b)
end

#%% Plot the flipping probability
# plot(omega_array, prob_array1, "-", label="ω = 0.1 GHz")
# plot!(omega_array, prob_array2, "r--", label="ω = 0.01 GHz")

plot(omega_array, prob_array3, label="ω1 = 0.5 GHz", lw=1.5)
plot!(omega_array, prob_array1, label="ω1 = 0.1 GHz", lw=1.5)
plot!(omega_array, prob_array2, label="ω1 = 0.01 GHz", lw=1.5)
xlabel!("ω (GHz)")
ylabel!("Flipping probability")

#%%
#####################################
# Fourier transform to observe inherent frequencies in probability curve

function power_spectrum(u0, tspan, p)
    p_b = prob_orientation(bloch!, u0, tspan, p)[2]
    t_array = prob_orientation(bloch!, u0, tspan, p)[1]
    fft_pb = real.(rfft(p_b))
    power_spec = broadcast(abs, fft_pb).^2
    power_spec[1] = 0
    n = size(t_array)[1]
    frequencies = rfftfreq(n, 0.1)
    omegaPrime = 2*π*frequencies
    return omegaPrime, power_spec
end

#%%
aR = 1.0; aL = 0.0; bR = 0.0; bL = 0.0
u0 = [aR,aL,bR,bL]

#%%
ω0 = 1.0; ω1 = 0.10
tmax = 10000.0
tspan = (0.0,tmax)

#%%
ω = 0.1
p = [ω, ω0, ω1]
omegaPrime1, Y1 = eval_power_spectrum(u0, tspan, p)

plt.title("Fourier transform of probability p_b")
plot(omegaPrime1, Y1, xlims = (0.0, 0.02))
xlabel!("ω' (GHz)")
ylabel!("Power spectrum")
plt.xlim([0, 1.4])

#%%
ω = 0.2
p = [ω, ω0, ω1]
omegaPrime2, Y2 = eval_power_spectrum(u0, tspan, p)


#%%
plot(omegaPrime2, Y2)

#%%
omegaPrime3, Y3 = eval_power_spectrum(u0, tspan, [0.6, 1.0, 1.0])
omegaPrime4, Y4 = eval_power_spectrum(u0, tspan, [0.7, 1.0, 1.0])
omegaPrime5, Y5 = eval_power_spectrum(u0, tspan, [0.8, 1.0, 1.0])
omegaPrime6, Y6 = eval_power_spectrum(u0, tspan, [0.9, 1.0, 1.0])

plot(omegaPrime3, Y3)
plot(omegaPrime4, Y4)
plot(omegaPrime5, Y5)
plot(omegaPrime6, Y6)
