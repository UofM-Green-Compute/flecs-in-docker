import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 
import os

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/SPH_Dust.txt")

# Define number of data points variable, arrays to store components

# Read data in from file
for i,row in enumerate(file): 
    if (i == 0):
        row = row.split("|") # Seperate no steps and no particles
        no_steps = int(row[0])
        no_particles = int(row[1])

        x_matrix = np.zeros((no_steps+1,no_particles))
        y_matrix = np.zeros((no_steps+1,no_particles))
        x_vel_matrix = np.zeros((no_steps+1,no_particles))
        y_vel_matrix = np.zeros((no_steps+1,no_particles))

    elif(i == 1):
        print(row)
    else:
        row = row.split("|")
        for j,pos_vec in enumerate(row):
            pos_vector = row[j].split(",")
            
            if(pos_vector[0] != "\n" and pos_vector[0] != ""):
                x_matrix[i-2][j] = float(pos_vector[0])
                y_matrix[i-2][j] = float(pos_vector[1])
                x_vel_matrix[i-2][j] = float(pos_vector[2])
                y_vel_matrix[i-2][j] = float(pos_vector[3])

# Graph formatting
fig, ax_ani = plt.subplots()
fig2, ax_hist = plt.subplots()
plt.suptitle("Coupled Oscillator")

ax_ani.set_title("Particles")
ax_ani.set_xlim(np.min(x_matrix)-1,np.max(x_matrix))
ax_ani.set_ylim(np.min(y_matrix),np.max(y_matrix))

fig.tight_layout()

# Animation
graph_particles, = ax_ani.plot([],[],'o',color="Red",label = 'Ideal')
_, _, graph_hist, = ax_hist.hist([],bins=50,density=True)

def particle(frame):

    state_x = x_matrix[frame,:]
    state_y = y_matrix[frame,:]

    graph_particles.set_data(state_x, state_y)

    return graph_particles

def histogram(frame):

    state_vx = x_vel_matrix[frame,:]
    state_vy = y_vel_matrix[frame,:]
    speed = np.sqrt( state_vx**2 + state_vy**2 )

    ax_hist.cla()
    _, _, graph_hist, = ax_hist.hist(speed,bins=50,fc="green",density=True)

    return graph_hist

ani = FuncAnimation(fig, particle, frames=no_steps, interval=1)
#ani2 = FuncAnimation(fig2, histogram, frames=no_steps, interval=1)

plt.show()