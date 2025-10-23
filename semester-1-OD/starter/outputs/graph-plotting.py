import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/starter/outputs/SHM-Data.txt")

# Define number of data points variable, arrays to store components
no_points = 0
pos = []
vel = []
acc = []

# Read data in from file
for i,row in enumerate(file): 
    if (i == 0):
        no_points = int(row)
    else:
        row = row.split(",")
        pos.append(float(row[0]))
        vel.append(float(row[1]))
        acc.append(float(row[2])) 

t = np.linspace(0,no_points,no_points)
zeroes = np.zeros((no_points,1))

# Graph formatting
fig, [ax, ax1, ax2, ax3] = plt.subplots(4)
plt.suptitle("Simple Harmonic Oscillator")
ax.set_title("Position")
ax1.set_title("Velocity")
ax2.set_title("Acceleration")
ax3.set_title("Particle")
fig.tight_layout()

# Animation
graph0, = ax.plot(t,pos,label="Position")
graph1, = ax1.plot(t,vel,color="Orange",label="Velocity")
graph2, = ax2.plot(t,acc,color="Red",label="Velocity")
graph3, = ax3.plot(0,0,'o',color="Red",label="Velocity")

ax3.set_xlim(min(pos),max(pos))
ax3.set_ylim(min(pos),max(pos))

for ax in fig.get_axes():
   ax.label_outer()        

def position(frame):
        return graph0.set_data(t[:frame],pos[:frame])

def velocity(frame):
    return graph1.set_data(t[:frame],vel[:frame])

def acceleration(frame):
    return graph2.set_data(t[:frame],acc[:frame])

def particle(frame):
    return graph3.set_data([pos[frame]],zeroes[frame])

ani = FuncAnimation(fig, position, frames=len(t), interval=10)
ani1 = FuncAnimation(fig, velocity, frames=len(t), interval=10)
ani2 = FuncAnimation(fig, acceleration, frames=len(t), interval=10)
ani3 = FuncAnimation(fig, particle, frames=len(t), interval=10)

plt.show()

# Making a gif
metadata = dict(title="Movie", artist="Oluwole")
writer = PillowWriter(fps=100, metadata=metadata)

time_list = []
pos_list = [] 
vel_list = []
acc_list = []
zero = []

with writer.saving(fig, "Oscillator.gif", 100):
    for i,time in enumerate(t):
        time_list.append(time)
        pos_list.append(pos[i])
        vel_list.append(vel[i])
        acc_list.append(acc[i])
        zero.append(0)

        graph0.set_data(time_list,pos_list)
        graph1.set_data(time_list,vel_list)
        graph2.set_data(time_list,acc_list)
        graph3.set_data(pos_list[-1:],zero[-1:])
        writer.grab_frame()


