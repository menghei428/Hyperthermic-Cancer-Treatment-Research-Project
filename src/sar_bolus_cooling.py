# 1D Pennes bioheat model for superficial hyperthermia
# FTCS in time + flux-conservative diffusion in space


# Referencing Literature
# https://pmc.ncbi.nlm.nih.gov/articles/PMC12835268/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC12835268/table/Tab2/


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

# Global settings

D = 0.04  # Total tissues depth in meters
treatment_time = 3600  # second

Nz = 201 # Number of spatial grid points
z = np.linspace(0, D, Nz)
dz = z[1] - z[0] # Spatial step size

dt = 0.02 # Time step (seconds)
Nt = int(treatment_time / dt)


# Parameters

rho_b = 1060 # blood density
c_b = 3770 # blood specific heat
T_body = 37

k_arr = np.zeros(Nz)         # thermal conductivity [W/(m K)]
rho_arr = np.zeros(Nz)       # density [kg/m^3]
cp_arr = np.zeros(Nz)        # specific heat [J/(kg K)]
omega_arr = np.zeros(Nz)     # blood perfusion [1/s]
qm_arr = np.zeros(Nz)        # metabolic heat [W/m^3]
Qapp_arr = np.zeros(Nz)      # applicator heating term [W/m^3]


# Layer

skin_end = 0.002    # 2 mm
fat_end = 0.015     # 15 mm
# muscle: remainder


# Water bolus cooling

h = 100  # W/(m^2 K),  heat transfer coefficient between bolus and skin. The literature suggested 70-152
T_bolus = 40 # bolus temperature


# Therapeutic threshold

threshold = 40

# Tissue properties

for i in range(Nz):
    zi = z[i]

    if zi < skin_end:  # skin
        k_arr[i] = 0.42
        rho_arr[i] = 1109
        cp_arr[i] = 3500
        omega_arr[i] = 0.0022
        qm_arr[i] = 1620
        Qapp_arr[i] = 1e4

    elif zi < fat_end:  # fat
        k_arr[i] = 0.25
        rho_arr[i] = 911
        cp_arr[i] = 2500
        omega_arr[i] = 0.00045
        qm_arr[i] = 300
        Qapp_arr[i] = 2e4

    else:  # muscle
        k_arr[i] = 0.5
        rho_arr[i] = 1090.4
        cp_arr[i] = 3600
        omega_arr[i] = 0.001
        qm_arr[i] = 480
        Qapp_arr[i] = 2e4

rhoCp = rho_arr * cp_arr
beta = (omega_arr * rho_b * c_b) / rhoCp
metabolic = qm_arr / rhoCp
Qapp_norm = Qapp_arr / rhoCp


#Stability check

alpha = k_arr / rhoCp
stability = np.max(alpha) * dt / dz**2

print(f"Stability: {stability:.2f}")

if stability > 0.5:
    print("Error, reduce dt or increase dz")
    sys.exit()


# Initial condition

T = np.ones(Nz) * T_body # Initial temperature everywhere = body temperature
T_new = T.copy()


# Harmonic mean conductivity at cell faces

# k_face[j] sits between node j and j+1
k_face = 2.0 * k_arr[:-1] * k_arr[1:] / (k_arr[:-1] + k_arr[1:])


# Time stepping

for n in range(Nt):

    # Interior nodes: flux-conservative diffusion

    flux_r = k_face[1:] * (T[2:] - T[1:-1]) / dz
    flux_l = k_face[:-1] * (T[1:-1] - T[:-2]) / dz

    diffusion = dt * (flux_r - flux_l) / (rhoCp[1:-1] * dz)
    perfusion = -beta[1:-1] * dt * (T[1:-1] - T_body)

    T_new[1:-1] = T[1:-1] + diffusion + perfusion + metabolic[1:-1] * dt + Qapp_norm[1:-1] * dt

    # Surface boundary at z = 0:
    # -k dT/dz = h (T_surface - T_bolus)
    # Use ghost node with central difference
    k0 = k_arr[0]
    ghost = T[1] - 2.0 * dz * (h / k0) * (T[0] - T_bolus)

    diff0 = k0 * dt * (T[1] - 2.0 * T[0] + ghost) / (rhoCp[0] * dz**2)
    perf0 = -beta[0] * dt * (T[0] - T_body)

    T_new[0] = (T[0] + diff0 + perf0 + metabolic[0] * dt + Qapp_norm[0] * dt)

    T_new[-1] = T_body

    T[:] = T_new[:]


# Find deepest point where temperature exceeds threshold

z_mm = z * 1000

if np.any(T >= threshold):
    treatable_depth = z_mm[T >= threshold][-1]
else:
    treatable_depth = 0

print(f"Final treatable depth (T >= {threshold:.1f} °C): {treatable_depth:.2f} mm")
print(f"Surface temperature: {T[0]:.2f} °C")
print(f"Max temperature: {np.max(T):.2f} °C at depth {z_mm[np.argmax(T)]:.2f} mm")

# Plotting graph

plt.figure(figsize=(8, 5))
plt.plot(z_mm, T, linewidth=2, color='darkred', label='Tissue temperature')
plt.axhline(threshold, linestyle='--', color='blue', label=f'{threshold:.0f} °C threshold')

if treatable_depth > 0:
    plt.axvline(treatable_depth, linestyle=':', color='purple',
                label=f'Treatable depth = {treatable_depth:.2f} mm')

# Layer shading
plt.axvspan(0, skin_end * 1000, color='pink', alpha=0.25, label='Skin')
plt.axvspan(skin_end * 1000, fat_end * 1000, color='orange', alpha=0.20, label='Fat')
plt.axvspan(fat_end * 1000, D * 1000, color='lightblue', alpha=0.15, label='Muscle')

plt.fill_between(z_mm, threshold, T, where=(T >= threshold),
                 color='green', alpha=0.20, label='Therapeutic zone')

plt.xlabel('Depth (mm)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature distribution after 1 hour of hyperthermia')
plt.xlim(0, D * 1000)
plt.ylim(T_body, max(np.max(T) + 0.5, T_bolus + 0.5))
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()