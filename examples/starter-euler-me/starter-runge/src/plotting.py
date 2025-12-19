import matplotlib.pyplot as plt
import numpy as np
import os

# Open data file
root_folder = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(root_folder, "outputs", "Coupled_Oscillators.txt")
save_path1 = os.path.join(root_folder, "outputs", "Coupled_Oscillators.pdf")
save_path2 = os.path.join(root_folder, "outputs", "Coupled_Oscillators.png")

file = open(data_path)

# Define number of data points variable, arrays to store components
no_points = 0
time = []

pos1 = []
vel1 = []
acc1 = []

pos2 = []
vel2 = []
acc2 = []

positions = [pos1,pos2]
velocities = [vel1,vel2]
accelerations = [acc1,acc2]

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
        pos1.append(float(row[1]))
        vel1.append(float(row[2]))
        acc1.append(float(row[3])) 

        pos2.append(float(row[4]))
        vel2.append(float(row[5]))
        acc2.append(float(row[6])) 

# Graph formatting
fig, axs = plt.subplots(3)

axs[0].set_xlim(0,float(max(time)))
axs[0].set_ylim(min(min(pos1),min(pos2))-1,max(max(pos1),max(pos2))+1)
axs[0].set_xlabel('time (s)')
axs[0].set_ylabel("Position (cm)")

axs[1].set_xlim(0,float(max(time)))
axs[1].set_ylim(min(min(vel1),min(vel2))-1,max(max(vel1),max(vel2))+1)
axs[1].set_ylabel("Velocity (cm s-1)")

axs[2].set_xlim(0,float(max(time)))
axs[2].set_ylim(min(min(acc1),min(acc2))-1,max(max(acc1),max(acc2))+1)
axs[2].set_ylabel("Acceleration (cm s-2)")

handle_1 = axs[0].plot(time, pos1, color = "tab:red", label="Particle 1")
axs[1].plot(time, vel1, color = "tab:red")
axs[2].plot(time, acc1, color = "tab:red")

handle_2 = axs[0].plot(time, pos2, color = "k", label="Particle 2")
axs[1].plot(time, vel2, color = "k")
axs[2].plot(time, acc2, color = "k")

for ax in fig.get_axes():
   ax.label_outer()  

handles = handle_1 + handle_2
labels = [handle.get_label() for handle in handles]
#fig.legend(handles, labels, ncol = 2, loc='lower center')
fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.05))
fig.align_ylabels(axs)
plt.tight_layout()
fig.savefig(save_path1, dpi=300, bbox_inches='tight')
fig.savefig(save_path2, dpi=300, bbox_inches='tight')
plt.show()
