from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# Material Variables
mass = 3
k = 4 * np.pi**2
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

    time_analytical = np.linspace(0, np.max(time), 10000)

    x1_analytical = x1model(time_analytical)
    x2_analytical = x2model(time_analytical)

    v1_analytical = v1model(time_analytical)
    v2_analytical = v2model(time_analytical)

    a1_analytical = a1model(time_analytical)
    a2_analytical = a2model(time_analytical)

    # Graph formatting
    # Formatting preamble
    cm = 1/2.54 #cm in inches
    mpl.rcParams['text.usetex'] = True
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
