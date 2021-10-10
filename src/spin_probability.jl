# Define the coupled ODEs
function bloch!(du, u, p, t)
    du[1] = 0.5 * (p[2]*u[2] + p[3] * cos(p[1]*t) * u[4])
    du[2] = -0.5 * (p[2] * u[1] + p[3] * cos(p[1]*t) * u[3])
    du[3] = -0.5 * (p[2] * u[4] - p[3] * cos(p[1]*t) * u[2])
    du[4] = 0.5 * (p[2] * u[3] - p[3] * cos(p[1]*t) * u[1])
end

function spin_down_probability(func, u0, tspan, p)
    problem = ODEProblem(func,u0,tspan,p)
    sol = solve(problem,RK4())
    prob = sol[3,:].^2 + sol[4,:].^2
    return sol.t, prob
end

#%%
## Plotting flipping probability
N = 1000
ω0 = 1.0; ω1 = 0.1
u0 = [1.0, 0, 0, 0]
tmax = 2.1 * pi / ω1
tspan = (0.0,tmax)
omega_array = range(0.9, 1.1, length=N)
prob_array1 = zeros(N)

for i = 1:N
    ω = omega_array[i]
    p = [ω, ω0, ω1]
    p_b = spin_down_probability(bloch!, u0, tspan, p)[2]
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
    p_b = spin_down_probability(bloch!, u0, tspan, param)[2]
    prob_array2[i] = maximum(p_b)
end

#%% Plot the flipping probability
plt.plot(omega_array, prob_array1, "-", label="ω = 0.1 GHz")
plt.plot(omega_array, prob_array2, "r--", label="ω = 0.01 GHz")
plt.xlabel("ω (GHz)")
plt.ylabel("Flipping probability")
plt.legend()

#%%
#####################################
using FFTW

function eval_power_spectrum(u0, tspan, param)
    p_b = spin_down_probability(bloch!, u0, tspan, param)[2]
    t_array = spin_down_probability(bloch!, u0, tspan, param)[1]
    fft_pb = fft(p_b)
    power_spec = broadcast(abs, fft_pb).^2
    power_spec[1] = 0
    n = size(t_array)[1]
    frequencies = fftfreq(n)
    omegaPrime = 2*pi*frequencies
    return omegaPrime, power_spec
end

#%%
aR = 1.0; aL = 0.0; bR = 0.0; bL = 0.0
u0 = [aR;aL;bR;bL]

#%%
ω0 = 1.0; ω1 = 0.10
tmax = 10000
tspan = (0.0,tmax)

#%%
ω = 0.1
param = [ω; ω0; ω1]
omegaPrime1, Y1 = eval_power_spectrum(u0, tspan, param)

plt.plot(omegaPrime1, Y1)
plt.title("Fourier transform of probability p_b")
plt.xlabel("ω (GHz)")
plt.ylabel("Power spectrum")
plt.xlim([0, 1.4])

#%%
ω = 0.2
param = [ω; ω0; ω1]
omegaPrime2, Y2 = eval_power_spectrum(u0, tspan, param)


#%%
plt.plot(omegaPrime2, Y2)
