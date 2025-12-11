import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 
import os
import scipy.constants as const

"""
Graph plotting for smoothed particle hydrodynamics code

Some potential additions:
Contouring irregularly spaced data: 
    - https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html
"""

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/SPH_Dust.txt")
specs_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/sph_code_specifications.txt"); 
density_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/density_field.txt")
density_gambel = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/gambel_density_field.txt")

# Read data in from file
for i,row in enumerate(specs_file): 
    if (i == 0):
        row = row.split("|") 
        no_steps = int(row[0])
        no_particles = int(row[1])
    if (i == 1):
        row = row.split("|") 
        Temp = int(row[0])
        Mass = float(row[1])
    if (i == 2):
        row = row.split("|") 
        GX = int(row[0])
        GX2 = int(row[1])
        GY = int(row[2])

# Define matrices to store positions
x_matrix = np.zeros((no_steps,no_particles))
y_matrix = np.zeros((no_steps,no_particles))
x_vel_matrix = np.zeros((no_steps,no_particles))
y_vel_matrix = np.zeros((no_steps,no_particles))
density_matrix = np.zeros((GY+1,(GX+GX2)+1))
gambel_matrix = np.zeros((15,2))

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

for i,row in enumerate(density_file): 
    row = row.split(",")
    for j,rho in enumerate(row):
        if(i < 5):
            density_matrix[i][j] = float(rho)

for i,row in enumerate(density_gambel): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix[i][j] = float(rho)

print(gambel_matrix)

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

# Graph formatting
fig, ax_ani = plt.subplots()
fig2, ax_hist = plt.subplots()
fig3, ax_density = plt.subplots()
fig4, ax_gambel = plt.subplots()
ax_ani.set_title("Particles")

ax_gambel.scatter(gambel_matrix[:,0], gambel_matrix[:,1])

# Contour plot for density 
density = ax_density.contourf(density_matrix)
fig3.colorbar(density,ax=ax_density)

# Edit number of bins for histogram
no_hist_bins = int(no_particles / 3) 

# Draw box
ax_ani.set_xlim(np.min(x_matrix)-1,np.max(x_matrix)+1)
ax_ani.set_ylim(np.min(y_matrix)-1,np.max(y_matrix)+1)
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

fig.tight_layout()

# Animation
graph_particles, = ax_ani.plot([],[],'o',color="Red",label = 'Ideal gas')
graph_maxwell, = ax_hist.plot(velocities, 2*10**6 * maxwell_boltzman_velocity(velocities),color = 'blue') # Scale of this plot is much smaller than hist
_, _, graph_hist_x, = ax_hist.hist([],bins=no_hist_bins,density=True,color='b')
ax_hist.scatter(v_rms, 2*10**6 * maxwell_boltzman_velocity(v_rms),color='red',marker='x')
# ax_hist.vlines(x = v_rms, ymin=0, ymax= 1*10**6 * max(maxwell_boltzman_velocity(velocities)) + 0.003, color = 'black', label = 'axvline - full height')

def particle(frame):
    state_x = x_matrix[frame,:]
    state_y = y_matrix[frame,:]

    graph_particles.set_data(state_x, state_y)

    return graph_particles

def histogram(frame):
    state_vx = x_vel_matrix[frame,:]
    state_vy = y_vel_matrix[frame,:]
    speed = np.sqrt( state_vx**2 + state_vy**2 )

    _, _, graph_hist_x, = ax_hist.hist(state_vx,bins=no_hist_bins,density=True,color='b')
    graph_maxwell, = ax_hist.plot(velocities, maxwell_boltzman_velocity(velocities))

ani = FuncAnimation(fig, particle, frames=no_steps, interval=0.001)
ani2 = FuncAnimation(fig2, histogram, frames=no_steps, interval=0.001)

plt.show()