# parameter sweep

import numpy as np
import matplotlib.pyplot as plt

def sensitivity_analysis(T_surface):

    #Properties

    D = 0.04  # Total tissues depth in meters
    rho_b = 1060 # blood density
    c_b = 3770 # blood specific heat
    T_body = 37
    treatment_time = 3600  # second

    # Slices

    Nz = 100 # Number of spatial grid points
    dz = D / (Nz - 1) # Spatial step size
    dt = 0.05 # Time step (seconds)
    Nt = int(treatment_time / dt) # Number of time steps: 72000

    #Grid

    z = np.linspace(0, D, Nz)

    # Parameter arrays

    alpha = np.zeros(Nz)
    beta = np.zeros(Nz)
    metabolic = np.zeros(Nz)



    for i in range(Nz):
        if z[i] < 0.002:  # skin
            k = 0.42 # thermal conductivity
            rho = 1109 # density
            C_p = 3500 # tissue specific heat
            omega_b = 0.0022 # blood perfusion rate
            q_m = 1620 # metabolic heat

        elif z[i] < 0.015:  # fat layer
            k = 0.25
            rho = 911
            C_p = 2500
            omega_b = 0.00045
            q_m = 300

        else:  # muscle
            k = 0.5
            rho = 1090.4
            C_p = 3600
            omega_b = 0.001
            q_m = 480

        alpha[i] = k / (rho * C_p)
        beta[i] = (omega_b * rho_b * c_b) / (rho * C_p)
        metabolic[i] = q_m / (C_p * rho)

    #Stability check

    stability = np.max(alpha) * dt / dz**2

    if stability > 0.5:
        raise ValueError(f"Unstable: CFL = {stability:.3f}. Reduce dt or increase dz.")


    # FTCS

    T = np.ones(Nz) * T_body # Initial temperature everywhere = body temperature
    T_new = T.copy()


    for i in range(Nt): # i: time step index
        diffusion = alpha[1:-1] * (dt / dz**2) * (T[2:] - 2*T[1:-1] + T[:-2])
        perfusion = -beta[1:-1] * dt * (T[1:-1] - T_body)
        
        T_new[1:-1] = T[1:-1] + diffusion + perfusion + metabolic[1:-1] * dt

        # Endpoints (excluded boundaries because the formula requires neighbouring points)

        T_new[0] = T_surface # The heat holds the skin at 43
        T_new[-1] = T_body # Neumann BC: No heat flux out of deep boundary

        T[:] = T_new[:]


        #if i % 200 == 0: # every 200 time steps => 360 frames
        #   T_history.append(T.copy())



    # Therapeutic threshold temperature
    threshold = 40

    # Find deepest point where temperature exceeds threshold

    z_mm = z * 1000

    if np.any(T >= threshold):
        treatable_depth = z_mm[T >= threshold][-1]
    else:
        treatable_depth = 0

    return treatable_depth

surface_temp = np.linspace(41, 47, 13)

depth_results = [] # treatable depth for each heat temp

for i in surface_temp:

    depth = sensitivity_analysis(i)

    depth_results.append(depth)

plt.figure()

plt.plot(surface_temp, depth_results,
         marker='o', linewidth=2)

plt.xlabel("Surface Heating Temperature (°C)")
plt.ylabel("Treatable Depth (mm)")
plt.title("Sensitivity of Treatable Depth to Surface Heating")

plt.grid(True)

plt.show()