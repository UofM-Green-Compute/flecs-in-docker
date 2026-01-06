from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
import os

# Material Variables
mass = 3
k = 4 * np.pi**2
omega = np.sqrt((3*k)/mass)
step_array = np.array([0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001])
r2_array = np.array([])

def model(t):
    x = 1-np.sqrt(mass/(3*k))*np.sin(omega*t)
    return x

for i, time_step in enumerate(step_array):
    # Open data file
    root_folder = os.path.dirname(os.path.dirname(__file__))
    data_path = os.path.join(root_folder, "outputs", f"Coupled_Oscillators_step={time_step:.6f}.txt")
    save_path1 = os.path.join(root_folder, "outputs", "Coupled_Oscillators.pdf")
    save_path2 = os.path.join(root_folder, "outputs", "Coupled_Oscillators.png")
    
    
    data = np.genfromtxt(data_path, skip_header=2, delimiter=",")
    time = data[:,0]
    pos1 = data[:,1]
    vel1 = data[:,2]
    acc1 = data[:,3]
    pos2 = data[:,4]
    vel2 = data[:,5]
    acc2 = data[:,6]

    x_model = model(time)
    r2 = r2_score(x_model, pos1)
    r2_array = np.append(r2_array, r2)
    print(r2)

    # Graph formatting
    fig, axs = plt.subplots(3)

    axs[0].set_xlim(0,2)
    axs[0].set_ylim(0,3)
    axs[0].set_xlabel('time (s)')
    axs[0].set_ylabel("Position (m)")

    axs[1].set_xlim(0,np.max(time))
    axs[1].set_ylim(np.min(vel1),np.max(vel2))
    axs[1].set_ylabel("Velocity (m s-1)")

    axs[2].set_xlim(0,np.max(time))
    axs[2].set_ylim(np.min(acc1),np.max(acc1))
    axs[2].set_ylabel("Acceleration (m s-2)")

    handle_1 = axs[0].plot(time, pos1, color = "tab:red", label="Particle 1")
    axs[1].plot(time, vel1, color = "tab:red")
    axs[2].plot(time, acc1, color = "tab:red")

    handle_2 = axs[0].plot(time, pos2, color = "k", label="Particle 2")
    axs[1].plot(time, vel2, color = "k")
    axs[2].plot(time, acc2, color = "k")

    handle_3 = axs[0].plot(time, x_model, color = "tab:red", linestyle = "--", label = "particle 1 analytical")
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
    plt.close()

save_path = os.path.join(root_folder, "outputs", "r2runge.txt")
save_data = np.column_stack((step_array, r2_array))
np.savetxt(save_path, save_data, delimiter=",")
plt.figure()
plt.scatter(np.log(step_array), np.log(1-r2_array))
plt.show()