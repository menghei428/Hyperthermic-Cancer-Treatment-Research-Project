# Pennes bioheat equation + metabolic heat generation
# FTCS method

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as Animation
import sys

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

# Parameter arrays

alpha = np.zeros(Nz)
beta = np.zeros(Nz)
metabolic = np.zeros(Nz)

#https://pmc.ncbi.nlm.nih.gov/articles/PMC12835268/table/Tab2/

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

# 1. Initialize placeholders
line, = ax.plot(z_mm, T_history[0], color='darkred', lw=2, label="Tissue Temp")
v_line = ax.axvline(0, color='purple', linestyle=':', label="Current Treatable Depth")
# Create an empty list to track our shading collection
fill_collection = [ax.fill_between(z_mm, threshold, T_history[0], where=(T_history[0] >= threshold), color='green', alpha=0.25)]

# 2. Set static visuals (Labels, Layers, and Threshold)
ax.set_xlim(0, 40)
ax.set_ylim(37, 43)
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Temperature (°C)")
ax.axhline(threshold, color='blue', linestyle='--', alpha=0.5, label="40°C Threshold")

# Visual Tissue Layers
ax.axvspan(0, 2, color='pink', alpha=0.2, label="Skin")
ax.axvspan(2, 15, color='orange', alpha=0.15, label="Fat")
ax.axvspan(15, 40, color='lightblue', alpha=0.1, label="Muscle")
ax.legend(loc='upper right', fontsize='small')

# 3. The Brain: The Update Function
def update(frame):
    # Get the temperature data for this specific time step
    current_T = T_history[frame]
    
    # A. Update the main temperature line
    line.set_ydata(current_T)
    
    # B. Calculate and update the Treatable Depth [cite: 10]
    if np.any(current_T >= threshold):
        current_depth = z_mm[current_T >= threshold][-1]
        v_line.set_xdata([current_depth, current_depth])
        v_line.set_visible(True)
    else:
        current_depth = 0
        v_line.set_visible(False)
        
    # C. Update the Shaded Therapeutic Zone [cite: 9]
    # We remove the old shading and add a new one for this frame
    global fill_collection
    fill_collection[0].remove()
    fill_collection[0] = ax.fill_between(z_mm, threshold, current_T, 
                                        where=(current_T >= threshold), 
                                        color='green', alpha=0.25)

    # D. Update Title (Fixes the 0.0 minutes bug!)
    current_time_mins = (frame * 200 * dt) / 60
    ax.set_title(f"Hyperthermia Evolution: {current_time_mins:.1f} mins\n"
                 f"Current Depth: {current_depth:.2f} mm")
    
    return line, v_line, fill_collection[0]

# 4. Run the Animation
# We set blit=False to ensure the Title and Fill_Between update correctly
ani = Animation.FuncAnimation(fig, update, frames=len(T_history), 
                              interval=50, blit=False, repeat=False)

plt.tight_layout()
plt.show()