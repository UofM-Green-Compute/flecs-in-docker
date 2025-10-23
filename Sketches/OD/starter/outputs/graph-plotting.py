import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation 

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/starter/outputs/Coupled_Oscillators.txt")

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
        print("Time (s), Position 1 (cm), Velocity 1 (cm s-1), " \
        "Acceleration 1 (cm s-2), Position 2 (cm), Velocity 2 (cm s-1), Acceleration 2 (cm s-2)")
    else:
        row = row.split(",")
        time.append(row[0])
        pos0.append(float(row[1]))
        vel0.append(float(row[2]))
        acc0.append(float(row[3])) 

        pos2.append(float(row[4]))
        vel2.append(float(row[5]))
        acc2.append(float(row[6])) 

t = np.linspace(0,no_points,no_points)
zeroes = np.zeros((no_points,1))

# Graph formatting
fig, [ax_pos, ax_vel, ax_acc, ax_ani] = plt.subplots(4)
plt.suptitle("Coupled Oscillator")
ax_pos.set_title("Position")
ax_pos.set_ylim(min(min(pos0),min(pos2))-1,max(max(pos0),max(pos2))+1)
xticks = np.arange(float(min(time)), float(max(time)), 5)
ax_vel.set_xticks(xticks)
ax_vel.set_title("Velocity")
#ax_vel.set_ylim(min(min(vel0),min(vel2)),max(max(vel0),max(vel2)))
ax_acc.set_title("Acceleration")
ax_acc.set_ylim(min(min(acc0),min(acc2)),max(max(acc0),max(acc2)))
ax_ani.set_title("Particle")
ax_ani.set_xlim(min(pos0)+1,max(pos0)+1)
ax_ani.set_ylim(min(pos0),max(pos0))
fig.tight_layout()

# Animation
graph0, = ax_pos.plot([],[],label="Position")
graph1, = ax_vel.plot(time,vel0,color="Orange",label="Velocity")
graph2, = ax_acc.plot([],[],color="Red",label="Velocity")
graph3, = ax_ani.plot(0,0,'o',color="Red",label="Velocity")

for ax in fig.get_axes():
   ax.label_outer()  

def all_positions(frame):
    for p in positions:
        line = ax_pos.plot(time[:frame+1],p[:frame+1],color='blue',label="Position")  

def all_velocities(frame):
    for v in velocities:
        line = ax_vel.plot(t[:frame+1],v[:frame+1],color='orange',label="Position") 

def all_accelerations(frame):
    for a in acceleration:
        line = ax_acc.plot(t[:frame+1],a[:frame+1],color='red',label="Position")  

def position(frame):
        return graph0.set_data(t[:frame],pos0[:frame])

def velocity(frame):
    return graph1.set_data(time[:frame],vel0[:frame])

def acceleration(frame):
    return graph2.set_data(t[:frame],acc0[:frame])

def particle(frame):
    return graph3.set_data([pos0[frame]],zeroes[frame])

ani = FuncAnimation(fig, all_positions, frames=len(t), interval=10)
ani1 = FuncAnimation(fig, velocity, frames=len(time), interval=10)
ani2 = FuncAnimation(fig, all_accelerations, frames=len(t), interval=10)
ani3 = FuncAnimation(fig, particle, frames=len(t), interval=10)

#ani = FuncAnimation(fig2, position2, frames=len(t), interval=10)
#ani1 = FuncAnimation(fig2, velocity2, frames=len(t), interval=10)
#ani2 = FuncAnimation(fig2, acceleration2, frames=len(t), interval=10)
#ani3 = FuncAnimation(fig2, particle2, frames=len(t), interval=10)

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
        pos_list.append(pos0[i])
        vel_list.append(vel0[i])
        acc_list.append(acc0[i])
        zero.append(0)

        graph0.set_data(time_list,pos_list)
        graph1.set_data(time_list,vel_list)
        graph2.set_data(time_list,acc_list)
        graph3.set_data(pos_list[-1:],zero[-1:])
        writer.grab_frame()


