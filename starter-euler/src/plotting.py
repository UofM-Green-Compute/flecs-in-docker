from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# Material Variables
mass = 3
k = 4 * np.pi**2
length = 1 # natural length of spring
omega = np.sqrt((3*k)/mass)
step_array = np.array([0.01, 0.005, 0.001])

def x1model(t):
    x = 1-np.sqrt(mass/(3*k))*np.sin(omega*t)
    return x

def x2model(t):
    x = 2+np.sqrt(mass/(3*k))*np.sin(omega*t)
    return x

def v1model(t):
    v = -omega*np.sqrt(mass/(3*k))*np.cos(omega*t)
    return v

def v2model(t):
    v = omega*np.sqrt(mass/(3*k))*np.cos(omega*t)
    return v

def a1model(t):
    a = omega**2*np.sqrt(mass/(3*k))*np.sin(omega*t)
    return a

def a2model(t):
    a = -omega**2*np.sqrt(mass/(3*k))*np.sin(omega*t)
    return a

# *** Make a 2x2 grid of plots *** 

# Open files
root_folder = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
euler_folder = os.path.dirname(os.path.dirname(__file__))
runge_folder = os.path.join(root_folder, "starter-runge")
Euler_t005_path = os.path.join(euler_folder, "outputs", f"Coupled_Oscillators_step={0.05:.6f}.txt") #E1
Euler_t001_path = os.path.join(euler_folder, "outputs", f"Coupled_Oscillators_step={0.01:.6f}.txt") #E2
Runge_t005_path = os.path.join(runge_folder, "outputs", f"Coupled_Oscillators_step={0.05:.6f}.txt") #E3
Runge_t001_path = os.path.join(runge_folder, "outputs", f"Coupled_Oscillators_step={0.01:.6f}.txt") #E4

# Save data from files into an array
data_E1 = np.genfromtxt(Euler_t005_path, skip_header=2, delimiter=",") #
time_E1 = data_E1[:,0]
pos_E1 = data_E1[:,1]
pos2_E1 = data_E1[:,4]

data_E2 = np.genfromtxt(Euler_t001_path, skip_header=2, delimiter=",")
time_E2 = data_E2[:,0]
pos_E2 = data_E2[:,1]
pos2_E2 = data_E2[:,4]

data_R1 = np.genfromtxt(Runge_t005_path, skip_header=2, delimiter=",") #
time_R1 = data_R1[:,0]
pos_R1 = data_R1[:,1]
pos2_R1 = data_R1[:,4]

data_R2 = np.genfromtxt(Runge_t001_path, skip_header=2, delimiter=",")
time_R2 = data_R2[:,0]
pos_R2 = data_R2[:,1]
pos2_R2 = data_R2[:,4]

# *** Plot graphs ***

# Graph formatting
cm = 1/2.54 #cm in inches
mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
fig2, axs2 = plt.subplots(2, 2, figsize=(17*cm, 10*cm), layout="constrained")

fig2.supxlabel("Time t/s")
fig2.supylabel("Position x/s")
for ax in axs2.flat:
        ax.label_outer()
        ax.tick_params(axis='both',direction='in')

axs2[0,0].set_ylim(0,3)
axs2[1,0].set_ylim(0,3)
axs2[0,1].set_ylim(0,3)
axs2[1,1].set_ylim(0,3)

# Graph plotting 
time_E1_analytical = np.linspace(0, np.max(time_E1), 10000)
x1_E1_analytical = x1model(time_E1_analytical)
x2_E1_analytical = x2model(time_E1_analytical)

time_E2_analytical = np.linspace(0, np.max(time_E2), 10000)
x1_E2_analytical = x1model(time_E2_analytical)
x2_E2_analytical = x2model(time_E2_analytical)

time_R1_analytical = np.linspace(0, np.max(time_R1), 10000)
x1_R1_analytical = x1model(time_R1_analytical)
x2_R1_analytical = x2model(time_R1_analytical)

time_R2_analytical = np.linspace(0, np.max(time_R2), 10000)
x1_R2_analytical = x1model(time_R2_analytical)
x2_R2_analytical = x2model(time_R2_analytical)

axs2[0,0].plot(time_E1, pos_E1, color = "k", linestyle = "--", label="Particle 1 Simulated")
axs2[0,0].plot(time_E1_analytical, x1_E1_analytical, color = "k", label = "Particle 1 Analytical")
axs2[0,0].plot(time_E1, pos2_E1, color = "tab:red", linestyle = "--", label="Particle 2 Simulated")
axs2[0,0].plot(time_E1_analytical, x2_E1_analytical, color = "tab:red", label = "Particle 2 Analytical")

axs2[1,0].plot(time_E2, pos_E2, color = "k", linestyle = "--")
axs2[1,0].plot(time_E2_analytical, x1_E2_analytical, color = "k")
axs2[1,0].plot(time_E2, pos2_E2, color = "tab:red", linestyle = "--")
axs2[1,0].plot(time_E2_analytical, x2_E2_analytical, color = "tab:red")

axs2[0,1].plot(time_R1, pos_R1, color = "k", linestyle = "--")
axs2[0,1].plot(time_R1_analytical, x1_R1_analytical, color = "k")
axs2[0,1].plot(time_R1, pos2_R1, color = "tab:red", linestyle = "--")
axs2[0,1].plot(time_R1_analytical, x2_R1_analytical, color = "tab:red")

axs2[1,1].plot(time_R2, pos_R2, color = "k", linestyle = "--")
axs2[1,1].plot(time_R2_analytical, x1_R2_analytical, color = "k")
axs2[1,1].plot(time_R2, pos2_R2, color = "tab:red", linestyle = "--")
axs2[1,1].plot(time_R2_analytical, x2_R2_analytical, color = "tab:red")

fig2.legend(loc='upper center',bbox_to_anchor = (0.5, 1.20),ncol = 2,frameon=False)

save_path = os.path.join(root_folder, "outputs", f"Coupled_Oscillators_step_t005_t001.pdf")
#fig2.savefig(save_path, bbox_inches = 'tight')

plt.show()

Energy_list = np.array([])
timeStep_list = np.array([])
fig1, ax1 = plt.subplots(figsize=(17*cm, 10*cm), layout="constrained")

fig1.supxlabel("Time t/s")
fig1.supylabel("Position x/s")
ax1.label_outer()
ax1.tick_params(axis='both',direction='in')
for i, time_step in enumerate(step_array):
    # Open data file
    root_folder = os.path.dirname(os.path.dirname(__file__))
    data_path = os.path.join(root_folder, "outputs", 
                             f"Coupled_Oscillators_step={time_step:.6f}.txt")
    save_path1 = os.path.join(root_folder, "outputs", 
                              f"Coupled_Oscillators_step={time_step:.6f}.pdf")
    save_path2 = os.path.join(root_folder, "outputs", 
                              f"Coupled_Oscillators_step={time_step:.6f}.png")
    
    data = np.genfromtxt(data_path, skip_header=2, delimiter=",")
    time = data[:,0]
    pos1 = data[:,1]
    vel1 = data[:,2]
    acc1 = data[:,3]
    pos2 = data[:,4]
    vel2 = data[:,5]
    acc2 = data[:,6]

    Energy = (mass*(vel1**2))/2 + (mass*(vel1**2))/2 + k*((pos1-length)**2+
                                                          (pos2-pos1-length)**2+
                                                          (2*length-pos2)**2)/2

    ax1.plot(Energy_list, pos_E1, color = "k", linestyle = "--", label="Particle 1 Simulated")


fig2.legend(loc='upper center',bbox_to_anchor = (0.5, 1.20),ncol = 2,frameon=False)
save_path = os.path.join(root_folder, "outputs", f"Coupled_Oscillators_step_t005_t001.pdf")
#fig2.savefig(save_path, bbox_inches = 'tight')
plt.show()

"""
    # Graph formatting
    # Formatting preamble
    cm = 1/2.54 #cm in inches
    # mpl.rcParams['text.usetex'] = True
    mpl.rcParams.update({'font.size': 14})
    fig, axs = plt.subplots(3, 1, figsize=(17*cm, 20*cm), layout="constrained")

    axs[0].set_xlim(0,2)
    axs[0].set_ylim(0,3)
    axs[0].set_ylabel("Position (m)")

    axs[1].set_xlim(0,np.max(time))
    axs[1].set_ylim(np.min(vel1),np.max(vel2))
    axs[1].set_ylabel("Velocity (m s-1)")

    axs[2].set_xlim(0,np.max(time))
    axs[2].set_ylim(np.min(acc1),np.max(acc1))
    axs[2].set_ylabel("Acceleration (m s-2)")

    axs[0].plot(time, pos1, color = "k", linestyle = "--", label="Particle 1: Computational")
    axs[1].plot(time, vel1, color = "k", linestyle = "--")
    axs[2].plot(time, acc1, color = "k", linestyle = "--")

    axs[0].plot(time, pos2, color = "tab:red", linestyle = "--", label="Particle 2: Computational")
    axs[1].plot(time, vel2, color = "tab:red", linestyle = "--")
    axs[2].plot(time, acc2, color = "tab:red", linestyle = "--")

    axs[0].plot(time_analytical, x1_analytical, color = "k", label = "Particle 1: Analytical")
    axs[1].plot(time_analytical, v1_analytical, color = "k")
    axs[2].plot(time_analytical, a1_analytical, color = "k")

    axs[0].plot(time_analytical, x2_analytical, color = "tab:red", label = "Particle 2: Analytical")
    axs[1].plot(time_analytical, v2_analytical, color = "tab:red")
    axs[2].plot(time_analytical, a2_analytical, color = "tab:red")

    handles, labels = axs[0].get_legend_handles_labels()

    leg = fig.legend(handles, labels, loc='upper center',
                 bbox_to_anchor=(0.5, 0), ncol = 2)
    for ax in axs.flat:
        ax.set(xlabel=r'time/s')
    for ax in axs.flat:
        ax.label_outer()
    fig.savefig(save_path1, bbox_inches = 'tight')
    fig.savefig(save_path2, bbox_inches = 'tight')
    plt.show()
"""