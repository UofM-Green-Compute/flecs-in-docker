import numpy as np
import matplotlib.pyplot as plt
import os

L = 0.01 # Length of system
Nx = 11 # Number of nodes
Pe = 1   # Peclet Number

#Boundary Conditions
Start = 273.15 # Temperature at start
End = 274.15   # Temperature at end

# Analytical Model as derived in lab book
def analytical_model(x, Pe_num, start, end):
    y = start + (((np.exp((Pe_num*x)/L))-1)/(np.exp(Pe_num)-1)) * (end - start)
    return y

# Open data file and generate array
root_folder = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(root_folder, "outputs", "temperature.txt")
save_pathpdf = os.path.join(root_folder, "outputs", f"Pe={Pe},Nodes={Nx}.pdf")
save_pathpng = os.path.join(root_folder, "outputs", f"Pe={Pe},Nodes={Nx}.png")
data = np.genfromtxt(data_path, delimiter = ",", skip_header = 1)

# Creates data according to the analytical model
x_model = np.linspace(0, L, num = 1000)
y_model = analytical_model(x_model, Pe, Start, End)

# Creates data according to the computational model
x_computational = data[:,0]
y_computational = data[:,1]

# Compares the analytical and computation model in a plot
plt.figure()
plt.plot(x_model, y_model, color = "k", linestyle = "-", label = "Exact")
plt.plot(x_computational, y_computational, color = "k", linestyle = "--", label = "Calculated")
ax = plt.gca()
# Set position and bounds of axes
ax.spines['bottom'].set_position(('data', Start))
ax.spines['top'].set_position(('data', End))
ax.spines['left'].set_bounds(Start, End)
ax.spines['right'].set_bounds(Start, End)
# Set ticks in correct position
ax.xaxis.set_ticks_position('bottom')
ax.xaxis.set_label_position('bottom')
ax.yaxis.set_ticks([tick for tick in ax.get_yticks() if (tick >= Start and tick <= End)])
plt.xlim(min(x_model), max(x_model))
plt.xlabel('x (m)')
plt.ylabel('T (K)')
#ax.yaxis.set_label_coords(-0.1, 0.6)
plt.legend(loc='upper left', bbox_to_anchor=(0, 0.95))
plt.savefig(save_pathpdf, dpi = 1200)
plt.savefig(save_pathpng, dpi = 1200)
plt.show()
plt.close()
