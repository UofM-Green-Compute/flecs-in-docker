import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

file = open("SHM-Data.txt")

no_points = 0
pos = []
vel = []
acc = []

for i,row in enumerate(file): 
    if (i == 0):
        no_points = int(row)
    else:
        row = row.split(",")
        pos.append(float(row[0]))
        vel.append(float(row[1]))
        acc.append(float(row[2])) 

t = np.linspace(0,no_points,no_points)

fig, ax = plt.subplots(3)
plt.suptitle("Simple Harmonic Oscillator")

ax[0].plot(t,pos,label="Position")
ax[0].set_title("Position")
ax[1].plot(t,vel,color='Orange',label="Velocity")
ax[1].set_title("Velocity")
ax[2].plot(t,acc,color='Red',label="Acceleration")
ax[2].set_title("Acceleration")
fig.tight_layout()

for ax in fig.get_axes():
    ax.label_outer()

artists = []

"""
ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=400)

fig, ax = plt.subplots()

def animate(i):
    ax.clear()
    point = pos[i]
    ax.plot(point,0,color = 'red', marker = 'o')
    ax.set_xlim([-5, 100])
    ax.set_ylim([-5, 5])

ani = animation.FuncAnimation(fig, animate, frames=no_points,interval=5, repeat=False)
#plt.close()
"""
plt.show()