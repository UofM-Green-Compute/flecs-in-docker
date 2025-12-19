import numpy as np
import matplotlib.pyplot as plt
import os

L = 1 # Length of system
Nx = 501 # Number of nodes
Pe = 50  # Peclet Number
Start = 0 # Conserved Quantity at start
End = 1 # Conserved Quantity at end

# Analytical Model as derived in lab book
def analytical_model(x, Pe_num, start, end):
    y = start + (((np.exp((Pe_num*x)/L))-1)/(np.exp(Pe_num)-1)) * (end - start)
    return y

# Open data file and generate array
root_folder = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(root_folder, "outputs", "1D_dynamic_fluid.txt")
save_pathpdf1 = os.path.join(root_folder, "outputs", f"Position_Pe={Pe},Nodes={Nx}.pdf")
save_pathpng1 = os.path.join(root_folder, "outputs", f"Position_Pe={Pe},Nodes={Nx}.png")
save_pathpdf2 = os.path.join(root_folder, "outputs", f"Time_Pe={Pe},Nodes={Nx}.pdf")
save_pathpng2 = os.path.join(root_folder, "outputs", f"Time_Pe={Pe},Nodes={Nx}.png")
data = np.genfromtxt(data_path, delimiter = ",", skip_header = 1)
STEPS = (len(data) - 1)

# Creates data according to the analytical model
x_model = np.linspace(0, 1, num = 1000)
phi_model = analytical_model(x_model, Pe, Start, End)

# Compares the analytical and computation model in a plot
plt.figure()
for i in range(11):
    data_plot = data[int(i*STEPS//10)]
    t = data[int(i*STEPS//10),0]
    # Creates data according to the computational model
    x_values = np.zeros_like(data_plot[1:])
    for index in range(len(x_values)):
        x_values[index] = index/(Nx-1)
    x_computational = (L * x_values)
    phi_computational = data_plot[1:]
    plt.plot(x_computational, phi_computational, color = "k", alpha = (i+1)/11, label = f"t = {t}s")
plt.plot(x_model, phi_model, color = "tab:red", linestyle = "-", label = "Exact Steady State")
ax = plt.gca()
# Set position and bounds of axes
ax.spines['bottom'].set_position(('data', 0))
ax.spines['top'].set_position(('data', 1))
ax.spines['left'].set_bounds(0, End)
ax.spines['right'].set_bounds(0, End)
# Set ticks in correct position
ax.xaxis.set_ticks_position('bottom')
ax.xaxis.set_label_position('bottom')
ax.yaxis.set_ticks([tick for tick in ax.get_yticks() if (tick >= 0 and tick <= End)])
plt.xlim(min(x_model), max(x_model))
plt.xlabel('x (m)')
plt.ylabel(r'$\phi$')
#ax.yaxis.set_label_coords(-0.1, 0.6)
plt.legend(loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.3))
plt.tight_layout()
plt.savefig(save_pathpdf1, dpi = 1200)
plt.savefig(save_pathpng1, dpi = 1200)
plt.show()
plt.close()

COlUMNS = len(data[0])
INDEX = COlUMNS//2
X = (INDEX - 1)/(Nx-1)
phi_asymptote = analytical_model(0.5, Pe, Start, End)
plt.figure()
plt.title(rf'$\phi(t)$ at x = {X}')
plt.plot(data[:,0], data[:,INDEX], color = 'k', label = 'unsteady state')
plt.hlines(phi_asymptote, 0, 0.5, colors='k', linestyles='--', label = 'steady state')
plt.xlabel("t (s)")
plt.ylabel(r"$\phi$")
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(save_pathpdf2, dpi = 1200)
plt.savefig(save_pathpng2, dpi = 1200)
plt.show()
plt.close()
