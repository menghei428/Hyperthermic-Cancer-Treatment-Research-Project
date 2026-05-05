# Pennes bioheat equation + metabolic heat generation
# FTCS method

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from tissue_properties import get_tissue_properties


#Properties

D = 0.04  # Total tissues depth in meters
rho_b = 1060 # blood density
c_b = 3770 # blood specific heat
T_body = 37
T_surface = 43
treatment_time = 3600  # second

# Slices

Nz = 100 # Number of spatial grid points
dz = D / (Nz - 1) # Spatial step size
dt = 0.05 # Time step (seconds)
Nt = int(treatment_time / dt) # Number of time steps: 72000

#Grid

z = np.linspace(0, D, Nz)

k_arr, rho_arr, cp_arr, omega_arr, qm_arr = get_tissue_properties(z)
rhoCp     = rho_arr * cp_arr
alpha     = k_arr / rhoCp
beta      = (omega_arr * rho_b * c_b) / rhoCp
metabolic = qm_arr / rhoCp

#Stability check

stability = np.max(alpha) * dt / dz**2

print(f"Stability: {stability:.2f}")

if stability > 0.5:
    print("Error, reduce dt or increase dz")
    sys.exit()


# FTCS

T = np.ones(Nz) * T_body # Initial temperature everywhere = body temperature
T_new = T.copy()
T_history = [] # For animation


for i in range(Nt): # i: time step index
    diffusion = alpha[1:-1] * (dt / dz**2) * (T[2:] - 2*T[1:-1] + T[:-2])
    perfusion = -beta[1:-1] * dt * (T[1:-1] - T_body)
    
    T_new[1:-1] = T[1:-1] + diffusion + perfusion + metabolic[1:-1] * dt

    # Endpoints (excluded boundaries because the formula requires neighbouring points)

    T_new[0] = T_surface # The heat holds the skin at 43
    T_new[-1] = T_body

    T[:] = T_new[:]

    if i % 200 == 0:
        T_history.append(T.copy())

T_history.append(T.copy())

# Therapeutic threshold temperature
threshold = 40

# Find deepest point where temperature exceeds threshold

z_mm = z * 1000

if np.any(T >= threshold):
    treatable_depth = z_mm[T >= threshold][-1]
else:
    treatable_depth = 0

print(f"Final treatable depth (mm): {treatable_depth:.2f} mm")


# Plotting graph

plt.plot(z_mm, T, linewidth=2, color='darkred', label="Tissue Temperature")
plt.axhline(threshold, linestyle='--', color='blue', label="40°C Threshold")

if treatable_depth > 0:
    plt.axvline(treatable_depth, linestyle=':', color='purple', label=f"Treatable Depth: {treatable_depth:.2f} mm")

# Skin layer regions
plt.axvspan(0, 2, color='pink', alpha=0.3, label="Skin")
plt.axvspan(2, 15, color='orange', alpha=0.2, label="Fat layer")
plt.axvspan(15, 40, color='lightblue', alpha=0.15, label="Muscle")

plt.fill_between(z_mm, threshold, T, where=(T >= threshold), color='green', alpha=0.25, label="Therapeutic Zone")
plt.xlabel("Depth (mm)")
plt.ylabel("Temperature (°C)")
plt.title("Final Temperature Distribution After 1 Hour (43°C)")
plt.xlim(0, 40)
plt.ylim(T_body, T_surface)
plt.legend()
plt.grid(True)
plt.show()


# Animation

fig, ax = plt.subplots(figsize=(8, 4))
z_mm = z * 1000 


line, = ax.plot(z_mm, T_history[0], color='darkred', lw=2, label="Tissue Temp")
v_line = ax.axvline(0, color='purple', linestyle=':', label="Current Treatable Depth")
# Create an empty list to track our shading collection
fill_collection = [ax.fill_between(z_mm, threshold, T_history[0], where=(T_history[0] >= threshold), color='green', alpha=0.25)]
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=9, va='top')


ax.set_xlim(0, 40)
ax.set_ylim(37, 43)
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Temperature (°C)")
ax.axhline(threshold, color='blue', linestyle='--', alpha=0.5, label="40°C Threshold")

ax.axvspan(0, 2, color='pink', alpha=0.2, label="Skin")
ax.axvspan(2, 15, color='orange', alpha=0.15, label="Fat")
ax.axvspan(15, 40, color='lightblue', alpha=0.1, label="Muscle")
ax.legend(loc='upper right', fontsize='small')


def update(frame):
    current_T = T_history[frame]
    
    line.set_ydata(current_T)
    
    if np.any(current_T >= threshold):
        current_depth = z_mm[current_T >= threshold][-1]
        v_line.set_xdata([current_depth, current_depth])
        v_line.set_visible(True)
    else:
        current_depth = 0
        v_line.set_visible(False)
        
    global fill_collection
    fill_collection[0].remove()
    fill_collection[0] = ax.fill_between(z_mm, threshold, current_T, 
                                        where=(current_T >= threshold), 
                                        color='green', alpha=0.25)

    current_time_mins = (frame * 200 * dt) / 60
    ax.set_title(f"Hyperthermia Evolution: {current_time_mins:.1f} mins\n"
                 f"Current Depth: {current_depth:.2f} mm")
    time_text.set_text(f"Frame: {frame} / {len(T_history)-1}, Time: {current_time_mins:.1f} min")
    
    return line, v_line, fill_collection[0], time_text

ani = animation.FuncAnimation(fig, update, frames=len(T_history), 
                              interval=50, blit=False, repeat=False)

plt.tight_layout()
ani.save('hyperthermia_simulation.gif', writer='pillow', fps=20)
plt.show()