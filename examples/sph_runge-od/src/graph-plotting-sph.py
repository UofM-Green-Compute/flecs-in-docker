import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 
import matplotlib.tri as tri
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import griddata

"""
Graph plotting for smoothed particle hydrodynamics code

Some potential additions:
Contouring irregularly spaced data: 
    - https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html
One colorbar
    - https://www.geeksforgeeks.org/data-visualization/how-to-have-one-colorbar-for-all-subplots-in-matplotlib/
"""

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_Dust.txt")
specs_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_specifications.txt"); 
density_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/density_field.txt")
density_gambel_20 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=20.txt")
density_gambel_50 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=50.txt")
density_gambel_75 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=75.txt")
density_gambel_100 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=100.txt")
density_gambel_200 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=200.txt")

# *** READ IN DATA ***

# Read Specifications file - info about the constants and parameters in the code 
for i,row in enumerate(specs_file): 
    if (i == 0):
        row = row.split("|")
        DT = float(row[0].split(":")[1])
        no_steps = int(row[1].split(":")[1])
        no_particles = int(row[2].split(":")[1])
    if (i == 1):
        row = row.split("|") 
        Temp = int(row[0].split(":")[1])
        Mass = float(row[1].split(":")[1])
    if (i == 2):
        row = row.split("|") 
        GX = int(row[0].split(":")[1])
        GX2 = int(row[1].split(":")[1])
        GY = int(row[2].split(":")[1])
    if (i == 3):
        row = row.split(":") 
        no_gambel_points = int(row[1])

# Define matrices to store positions of particles 
x_matrix = np.zeros((no_steps,no_particles))
y_matrix = np.zeros((no_steps,no_particles))
density_on_particle_matrix = np.zeros((no_steps,no_particles))
x_vel_matrix = np.zeros((no_steps,no_particles))
y_vel_matrix = np.zeros((no_steps,no_particles))
density_matrix = np.zeros((no_steps*(GY+1),(GX+GX2)+1))
density = np.zeros((GY+1,(GX+GX2)+1))
gambel_matrix_20 = np.zeros((no_gambel_points,2))
gambel_matrix_50 = np.zeros((no_gambel_points,2))
gambel_matrix_75 = np.zeros((no_gambel_points,2))
gambel_matrix_100 = np.zeros((no_gambel_points,2))
gambel_matrix_200 = np.zeros((no_gambel_points,2))
energy = np.zeros((no_steps,2))

# Read data on particle positions into matrices defined above
for i,row in enumerate(file): 
    if(i == 0):
        print(row)
    else:
        row = row.split("|")
        for j,pos_vec in enumerate(row):
            pos_vector = row[j].split(",")
            
            if(pos_vector[0] != "\n" and pos_vector[0] != ""):
                x_matrix[i-1][j] = float(pos_vector[0])
                y_matrix[i-1][j] = float(pos_vector[1])
                x_vel_matrix[i-1][j] = float(pos_vector[2])
                y_vel_matrix[i-1][j] = float(pos_vector[3])
                density_on_particle_matrix[i-1][j] = float(pos_vector[3])

# Read in data on density distribution
for i,row in enumerate(density_file): 
    row = row.split(",")
    for j,rho in enumerate(row):
            density_matrix[i][j] = float(rho)

# Read in data on how estimation error varies with h
for i,row in enumerate(density_gambel_20): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_20[i][j] = float(rho)

for i,row in enumerate(density_gambel_50): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_50[i][j] = float(rho)

for i,row in enumerate(density_gambel_75): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_75[i][j] = float(rho)

for i,row in enumerate(density_gambel_100): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_100[i][j] = float(rho)

for i,row in enumerate(density_gambel_200): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_200[i][j] = float(rho)

# *** FUNCTIONS ***

# Get density distribution at a particular time step
def get_density(time_step):
    density_m = np.zeros((GY+1,(GX+GX2)+1))
    for i in range(GY+1):
        for j in range(GX+GX2+1): 
            density_m[i][j] = density_matrix[i + time_step * (GY+1)][j]
    return density_m

# Plot the Maxwell-Boltzmann Distribution
k = 1.380649 * 10**-23
v_rms = (( k * Temp) / Mass ) ** 0.5
velocities = np.linspace(-1000, 1000, 2000)

def maxwell_boltzman_velocity(v):
    N = (Mass / (2 * np.pi * k * Temp))**1.5 
    return N * np.exp( - (Mass * v**2 ) / (2 * k * Temp) )

def maxwell_boltzman_speed(v):
    N = 4 * np.pi * (Mass / (2 * np.pi * k * Temp))**1.5 
    return N * v**2 * np.exp( - (Mass * v**2 ) / (2 * k * Temp) )

# Fit data to a straight line
def linear_model(x, m, c):
    return m*x + c

# Linear fit of gambel density
#logH_Gambel20 = np.log(gambel_matrix_20[:,0])
#logGambel20 = np.log(gambel_matrix_20[:,1])
#logH_Gambel50 = np.log(gambel_matrix_50[:,0])
#logGambel50 = np.log(gambel_matrix_50[:,1])

#param20, param_cov20 = curve_fit(linear_model, logH_Gambel20, logGambel20)
#y = param20[0]*logGambel20 + param20[1]
#print("HIII ",param20[0],param20[1])

# Calculate energy of the simulation at each time step
for i in range(no_steps):
    energy_val = 0
    for j in range(no_particles):
        energy_val += 1/2 * Mass * ( x_vel_matrix[i][j]**2 + y_vel_matrix[i][j]**2 ) 
    energy[i][0] = i 
    energy[i][1] = energy_val

# Got this function from youtuber https://youtu.be/TZvt_fD8VP4
def interpolate_data(x, y, values, new_size, method='linear'):
    # Flatten coordinates
    points = np.column_stack((x.flatten(), y.flatten()))
    
    # Create regular grid
    x_new = np.linspace(x.min(), x.max(), new_size)
    y_new = np.linspace(y.min(), y.max(), new_size)
    X_new, Y_new = np.meshgrid(x_new, y_new)
    
    # Interpolate
    values_new = griddata(points, values.flatten(), 
                         (X_new, Y_new), method=method)
    
    return X_new, Y_new, values_new

# *** FORMAT GRAPHS ***

# Create figures and axes 
cm = 1/2.54 #cm in inches
fig, ax_ani = plt.subplots(figsize=(17*cm,10*cm)) 
fig2, ax_hist = plt.subplots(figsize=(17*cm,10*cm))
fig3, ax_density = plt.subplots(figsize=(17*cm,10*cm))
fig3b, ax_density3b = plt.subplots(2,2,figsize=(17*cm,10*cm))
fig4, ax_gambel = plt.subplots(figsize=(17*cm,10*cm))
fig5, ax_energy = plt.subplots(figsize=(17*cm,10*cm))

# Format graphs
mpl.rcParams.update({'font.size': 12})

# ax_density.set_aspect('equal',adjustable='box')
ax_ani.set_aspect('equal',adjustable='box')

# Major grid:
ax_gambel.grid(True, which='major', linestyle='-', linewidth=0.75, alpha=0.25)
ax_hist.grid(True, which='major', linestyle='-', linewidth=0.75, alpha=0.25)
# Minor ticks and grid:
ax_gambel.minorticks_on()
ax_gambel.grid(True, which='minor', linestyle='-', linewidth=0.25, alpha=0.15)
ax_gambel.set_axisbelow(True) # <-- Ensure grid is below data
# ax_gambel.tick_params(axis='both',direction='in')

# Set figure titles
ax_ani.set_title("Particles")
ax_ani.set_xlim(np.min(x_matrix)-1,np.max(x_matrix)+1)
ax_ani.set_ylim(np.min(y_matrix)-1,np.max(y_matrix)+1)

ax_density.set_title("Contour plot of particle density")
ax_density.set_xlabel("Box Width")
ax_density.set_ylabel("Box Height")

ax_energy.set_title("Calculated system energy over time")
ax_energy.set_xlabel("Time")
ax_energy.set_ylabel("Percentage change in system energy ")

ax_gambel.set_xlabel("h")
ax_gambel.set_ylabel("Average Error")

fig.tight_layout()

# *** PLOT DATA ***

# --- DENSITY FIELD ---
# Select grid points at a 4 time steps
step_1 = 0
step_2 = 0
step_3 = 0
step_4 = 0

x_00 = x_matrix[step_1,:]
y_00 = y_matrix[step_1,:]
z_00 = density_on_particle_matrix[step_1,:]

x_01 = x_matrix[step_2,:]
y_01 = y_matrix[step_2,:]
z_01 = density_on_particle_matrix[step_2,:]

x_10 = x_matrix[step_3,:]
y_10 = y_matrix[step_3,:]
z_10 = density_on_particle_matrix[step_3,:]

x_11 = x_matrix[step_4,:]
y_11 = y_matrix[step_4,:]
z_11 = density_on_particle_matrix[step_4,:]

# -- Tricontour interpolation -- 
ax_density3b[0,0].tricontour(x_00, y_00, z_00, levels=14, linewidths=0.5, colors='k')
cntr2_00 = ax_density3b[0,0].tricontourf(x_00, y_00, z_00, levels=14, cmap="RdBu_r")
ax_density3b[0,0].plot(x_00, y_00, 'ko', ms=3)

ax_density3b[0,1].tricontour(x_01, y_01, z_01, levels=14, linewidths=0.5, colors='k')
cntr2_01 = ax_density3b[0,1].tricontourf(x_01, y_01, z_01, levels=14, cmap="RdBu_r")
ax_density3b[0,1].plot(x_01, y_01, 'ko', ms=3)

ax_density3b[1,0].tricontour(x_10, y_10, z_10, levels=14, linewidths=0.5, colors='k')
cntr2_10 = ax_density3b[1,0].tricontourf(x_10, y_10, z_10, levels=14, cmap="RdBu_r")
ax_density3b[1,0].plot(x_10, y_10, 'ko', ms=3)

ax_density3b[1,1].tricontour(x_11, y_11, z_11, levels=14, linewidths=0.5, colors='k')
cntr2_11 = ax_density3b[1,1].tricontourf(x_11, y_11, z_11, levels=14, cmap="RdBu_r")
ax_density3b[1,1].plot(x_11, y_11, 'ko', ms=3)

fig3b.colorbar(cntr2_00, ax=ax_density3b[0,0])
fig3b.colorbar(cntr2_01, ax=ax_density3b[0,1])
fig3b.colorbar(cntr2_10, ax=ax_density3b[1,0])
fig3b.colorbar(cntr2_11, ax=ax_density3b[1,1])

"""
# -- Grid Interpolation -- 
ngridx = 100
ngridy = 200
xi = np.linspace(0, GX+GX2, ngridx)
yi = np.linspace(0, GY, ngridy)
triang = tri.Triangulation(x, y)
interpolator = tri.LinearTriInterpolator(triang, z)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)
ax_density3b[0,1].contour(xi, yi, zi, levels=14, linewidths=0.5, colors='k')
cntr1 = ax_density3b[0,1].contourf(xi, yi, zi, levels=14, cmap="RdBu_r")
fig3b.colorbar(cntr1, ax=ax_density3b)
# -- Youtuber interpolation -- 
x_int, y_int, z_int = interpolate_data(x,y,z,100)
ax_density3b[1,1].contourf(x_int,y_int,z_int)
"""

# Energy 
ax_energy.plot(energy[:,0],((energy[:,1] - energy[0][1])/energy[0][1])*100) # Percentage change in energy vs time
ax_energy.hlines(y = 0, xmin=min(energy[:,0]), xmax=max(energy[:,0]), color = 'black', linestyle='--')
E_KT = k * Temp 
# ax_energy.plot(energy[:,0],energy[:,1])
# ax_energy.hlines(y = E_KT , xmin=min(energy[:,0]), xmax=max(energy[:,0]), color = 'red', linestyle='--')

# Scatter plot of Gambel errors 
#ax_gambel.scatter(gambel_matrix_20[:,0], gambel_matrix_20[:,1], s=15, color='k', marker='.',label = '20 Samples')
ax_gambel.scatter(gambel_matrix_50[:,0], gambel_matrix_50[:,1], s=15, color='k', marker='^',label = '50 Samples')
# ax_gambel.scatter(gambel_matrix_75[:,0], gambel_matrix_75[:,1], s=15, color='k', marker=',',label = '75 Samples')
ax_gambel.scatter(gambel_matrix_100[:,0], gambel_matrix_100[:,1], s=15, color='k', marker=',',label = '100 Samples')
ax_gambel.scatter(gambel_matrix_200[:,0], gambel_matrix_200[:,1], color='k', marker='o',label = '200 Samples')
# ax_gambel.scatter(logH_Gambel20, logGambel20)
# ax_gambel.scatter(logH_Gambel50, logGambel50)

print(min(gambel_matrix_20[:,1]))
print(min(gambel_matrix_50[:,1]))
print(min(gambel_matrix_75[:,1]))
print(min(gambel_matrix_100[:,1]))
print(min(gambel_matrix_200[:,1]))

# Edit number of bins for histogram
no_hist_bins = int(no_particles / 3) 

# Draw box
ax_ani.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 5, ymin=0, ymax=2, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 5, ymin=3, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 0, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 5, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')

ax_density.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_density.vlines(x = 5, ymin=0, ymax=2, color = 'black', label = 'axvline - full height')
ax_density.vlines(x = 5, ymin=3, ymax=5, color = 'black', label = 'axvline - full height')
ax_density.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_density.hlines(y = 0, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')
ax_density.hlines(y = 5, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')

for ax in ax_density3b.flatten():
    ax.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
    ax.vlines(x = 5, ymin=0, ymax=2, color = 'black', label = 'axvline - full height')
    ax.vlines(x = 5, ymin=3, ymax=5, color = 'black', label = 'axvline - full height')
    ax.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
    ax.hlines(y = 0, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')
    ax.hlines(y = 5, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')

# Generate graph for SPH example
"""
ax_ani.set_title("") 
ax_ani.set_xlim(0,5)
ax_ani.set_ylim(0,5)
ax_ani.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 5, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.plot(x_matrix[0,:],y_matrix[0,:],'o',color="Red",label = 'Smoothed Particle Hydrodynamics')
ax_ani.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 0, xmin=0, xmax=5, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 5, xmin=0, xmax=5, color = 'black', label = 'axvline - full height')
ax_ani.tick_params(axis='both', direction='in')
time_step = 0
particle_no=1
draw_circle = plt.Circle((x_matrix[time_step,particle_no],y_matrix[time_step,particle_no]),1,color='black',linestyle="--",fill=False)
x = np.cos(np.pi / 4)
y = np.cos(np.pi / 4)
ax_ani.arrow(x_matrix[time_step,particle_no], y_matrix[time_step,particle_no], x, y, head_width=0.15, color='black', length_includes_head=True)
ax_ani.annotate("h",xy=(x_matrix[time_step,particle_no]+x/2-0.3,y_matrix[time_step,particle_no]+y/2-0.1))
ax_ani.add_artist(draw_circle)
"""

# fig4.legend()
fig4.legend(loc='upper center',bbox_to_anchor = (0.5, 1.0),ncol = 3,frameon=False)

# Save figures
#fig.savefig("outputs/SPH_Schematic.pdf") 
fig3.savefig("outputs/Density_Contour.pdf")
fig3b.savefig("outputs/Density_Irregular_Contour.pdf")
fig4.savefig("outputs/Spline_Kernel_Estimation_Error.pdf")
fig5.savefig("outputs/System_Energy_Variation.pdf")

# Animation
graph_particles, = ax_ani.plot([],[],'o',color="Red",label = 'Ideal gas')
graph_maxwell, = ax_hist.plot(velocities, 2*10**6 * maxwell_boltzman_velocity(velocities),color = 'blue') # Scale of this plot is much smaller than hist
_, _, graph_hist_x, = ax_hist.hist([],bins=no_hist_bins,density=True,color="#4f41e7",alpha=0.7) 
graph_contour = ax_density.contourf(density)
ax_hist.scatter(v_rms, 2*10**6 * maxwell_boltzman_velocity(v_rms),color='red',marker='x')

def particle(frame):
    state_x = x_matrix[frame,:]
    state_y = y_matrix[frame,:]

    graph_particles.set_data(state_x, state_y)

    return graph_particles

def histogram(frame):
    state_vx = x_vel_matrix[frame,:]
    state_vy = y_vel_matrix[frame,:]
    speed = np.sqrt( state_vx**2 + state_vy**2 )

    _, _, graph_hist_x, = ax_hist.hist(state_vx,bins=no_hist_bins,density=True,color="#595898",alpha=0.7)
    graph_maxwell = ax_hist.plot(velocities, maxwell_boltzman_velocity(velocities))

# Contour plot for density 
def contour(frame):

    global graph_contour

    density = get_density(frame)

    for coll in ax_density.collections:
        coll.remove()

    ax_density.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
    ax_density.vlines(x = 5, ymin=0, ymax=2, color = 'black', label = 'axvline - full height')
    ax_density.vlines(x = 5, ymin=3, ymax=5, color = 'black', label = 'axvline - full height')
    ax_density.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
    ax_density.hlines(y = 0, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')
    ax_density.hlines(y = 5, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')

    graph_contour = ax_density.contourf(density)

    return graph_contour

ani = FuncAnimation(fig, particle, frames=no_steps, interval=0.001)
ani2 = FuncAnimation(fig2, histogram, frames=no_steps, interval=0.001)
ani3 = FuncAnimation(fig3, contour, frames=no_steps, interval=0.001)

plt.show()