import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 

# Open data file
file = open("Coupled_Oscillators.txt")

# Define number of data points variable, arrays to store components
no_points = 0
time = []

pos0 = []
vel0 = []
acc0 = []

pos2 = []
vel2 = []
acc2 = []

positions = [pos0,pos2]
velocities = [vel0,vel2]
accelerations = [acc0,acc2]

# Read data in from file
for i,row in enumerate(file): 
    if (i == 0):
        row = row.split(": ")
        no_points = int(row[1])
        print(no_points)
    elif(i == 1):
        print(row)
    else:
        row = row.split(",")
        time.append(float(row[0]))
        pos0.append(float(row[1]))
        vel0.append(float(row[2]))
        acc0.append(float(row[3])) 

        pos2.append(float(row[4]))
        vel2.append(float(row[5]))
        acc2.append(float(row[6])) 

zeroes = np.zeros((no_points,1))

# Graph formatting
fig, [ax_pos, ax_vel, ax_acc, ax_ani] = plt.subplots(4)
plt.suptitle("Coupled Oscillator")

ax_pos.set_title("Position")
ax_pos.set_xlim(0,float(max(time)))
ax_pos.set_ylim(min(min(pos0),min(pos2))-1,max(max(pos0),max(pos2))+1)

ax_vel.set_title("Velocity")
ax_vel.set_xlim(0,float(max(time))+1)
ax_vel.set_ylim(min(min(vel0),min(vel2))-1,max(max(vel0),max(vel2))+1)

ax_acc.set_title("Acceleration")
ax_acc.set_xlim(0,float(max(time))+1)
ax_acc.set_ylim(min(min(acc0),min(acc2))-1,max(max(acc0),max(acc2))+1)

ax_ani.set_title("Particle")
ax_ani.set_xlim(min(min(pos0),min(pos2))-1,max(max(pos0),max(pos2))+1)
ax_ani.set_ylim(-1,1)

fig.tight_layout()

# Animation
#   Position graphs
graph0, = ax_pos.plot([],[],color = 'blue',label="Particle 1")
graph0b, = ax_pos.plot([],[],color = 'blue',label="Particle 2",linestyle='--')
#   Velocity graphs
graph1, = ax_vel.plot([],[],color="Orange")
graph1b, = ax_vel.plot([],[],color="Orange",linestyle='--')
#   Acceleration graphs
graph2, = ax_acc.plot([],[],color="Red")
graph2b, = ax_acc.plot([],[],color="Red",linestyle='--')
#   Anmation graphs
graph3, = ax_ani.plot(0,0,'o',color="Red",label = 'Particle 1')
graph3b, = ax_ani.plot(0,0,'s',color="Blue",label = 'Particle 1')

for ax in fig.get_axes():
   ax.label_outer()  

# Make 2 legends
fig.legend([graph0, graph0b], ['Particle 1', 'Partcile 2'],
          loc='upper left', frameon=False)
fig.legend([graph3, graph3b], ['Particle 1', 'Partcile 2'],
          loc='upper right', frameon=False)

def all_positions(frame):
    graph0.set_data(time[:frame],pos0[:frame])
    graph0b.set_data(time[:frame],pos2[:frame])  

def all_velocities(frame):
    graph1.set_data(time[:frame],vel0[:frame])
    graph1b.set_data(time[:frame],vel2[:frame]) 

def all_accelerations(frame):
    graph2.set_data(time[:frame],acc0[:frame])
    graph2b.set_data(time[:frame],acc2[:frame])   

def particle(frame):
    graph3.set_data([pos0[frame]],zeroes[frame])
    graph3b.set_data([pos2[frame]],zeroes[frame])

ani = FuncAnimation(fig, all_positions, frames=len(time), interval=10)
ani1 = FuncAnimation(fig, all_velocities, frames=len(time), interval=10)
ani2 = FuncAnimation(fig, all_accelerations, frames=len(time), interval=10)
ani3 = FuncAnimation(fig, particle, frames=len(time), interval=10)

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
    for i,time in enumerate(time):
        time_list.append(time)
        pos_list.append(pos0[i])
        vel_list.append(vel0[i])
        acc_list.append(acc0[i])
        zero.append(0)

        graph0.set_data(time_list,pos_list)
        graph1.set_data(time_list,vel_list)
        graph2.set_data(time_list,acc_list)
        graph3.set_data(pos_list[-1:],zero[-1:])
        writer.grab_frame()