import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 
import os
import scipy.constants as const

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/SPH_Dust.txt")
specs_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/sph_code_specifications.txt"); 

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

x_matrix = np.zeros((no_steps,no_particles))
y_matrix = np.zeros((no_steps,no_particles))
x_vel_matrix = np.zeros((no_steps,no_particles))
y_vel_matrix = np.zeros((no_steps,no_particles))

no_hist_bins = int(no_particles )

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

# Maxwell-Boltzmann Distribution
k = 1.380649 * 10**-23
peak = ((2 * k * Temp) / Mass ) ** 0.5

x = np.linspace(0, 25, 100) 
velocities = np.linspace(0, 4000, 1000)

def maxwell_boltzman(v):
    N = 4 * np.pi * (v**2) * (Mass / (2 * np.pi * k * Temp))**1.5 
    return N * np.exp( - (Mass * v**2 ) / (2 * k * Temp))

# Graph formatting
fig, ax_ani = plt.subplots()
fig2, ax_hist = plt.subplots()

ax_ani.set_title("Particles")

ax_ani.set_xlim(np.min(x_matrix)-1,np.max(x_matrix)+1)
ax_ani.set_ylim(np.min(y_matrix)-1,np.max(y_matrix)+1)
ax_ani.vlines(x = 0, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 5, ymin=0, ymax=2, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 5, ymin=3, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.vlines(x = 10, ymin=0, ymax=5, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 0, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')
ax_ani.hlines(y = 5, xmin=0, xmax=10, color = 'black', label = 'axvline - full height')

fig.tight_layout()

# Animation
graph_particles, = ax_ani.plot([],[],'o',color="Red",label = 'Ideal')
# graph_maxwell, = ax_hist.plot(velocities, maxwell_boltzman(velocities)) # Scale of this plot is much smaller than hist
_, _, graph_hist, = ax_hist.hist([],bins=no_hist_bins,density=True,color='b')
# ax_hist.vlines(x = peak, ymin=0, ymax=max(maxwell_boltzman(velocities))+0.0001, color = 'black', label = 'axvline - full height')

def particle(frame):
    state_x = x_matrix[frame,:]
    state_y = y_matrix[frame,:]

    graph_particles.set_data(state_x, state_y)

    return graph_particles

def histogram(frame):
    state_vx = x_vel_matrix[frame,:]
    state_vy = y_vel_matrix[frame,:]
    speed = np.sqrt( state_vx**2 + state_vy**2 )

    _, _, graph_hist, = ax_hist.hist(speed,bins=no_hist_bins,density=True,color='b')

    # return graph_hist

ani = FuncAnimation(fig, particle, frames=no_steps, interval=0.001)
ani2 = FuncAnimation(fig2, histogram, frames=no_steps, interval=0.001)

plt.show()