import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

L = 0.01 # Length of system
N = np.array([41, 21, 11]) # Number of nodes
Pe = 20  # Peclet Number

#Boundary Conditions
Start = 273.15 # Temperature at start
End = 274.15   # Temperature at end

# Analytical Model as derived in lab book
def analytical_model(x, Pe_num, start, end):
    y = start + (((np.exp((Pe_num*x)/L))-1)/(np.exp(Pe_num)-1)) * (end - start)
    return y

# Open data file and generate array
root_folder = os.path.dirname(os.path.dirname(__file__))
data_path1 = os.path.join(root_folder, "outputs", "temperature,N=41.txt")
data_path2 = os.path.join(root_folder, "outputs", "temperature,N=21.txt")
data_path3 = os.path.join(root_folder, "outputs", "temperature,N=11.txt")
save_path = os.path.join(root_folder, "outputs", f"Pe={Pe}.pdf")
data1 = np.genfromtxt(data_path1, delimiter = ",", skip_header = 1)
data2 = np.genfromtxt(data_path2, delimiter = ",", skip_header = 1)
data3 = np.genfromtxt(data_path3, delimiter = ",", skip_header = 1)

# Creates data according to the analytical model
x_model = np.linspace(0, L, num = 1000)
y_model = analytical_model(x_model, Pe, Start, End)

# Creates data according to the computational model
x_computational1 = data1[:,0]
T_computational1 = data1[:,1]
x_computational2 = data2[:,0]
T_computational2 = data2[:,1]
x_computational3 = data3[:,0]
T_computational3 = data3[:,1]

comparison1_T1 = T_computational1[0::4][1:-1]
comparison2_T1 = T_computational1[0::2][1:-1]
comparison_T2 = T_computational2[0::2][1:-1]
comparison_T3 = T_computational3[0::2][1:-1]

# After running once with Pe = 1 it was confirmed p = 2
"""p_list = np.log((comparison_T2-T_computational3[1:-1])/(comparison1_T1-
                                                         comparison_T2))/np.log(2)
p = np.mean(p_list)"""
p = 2

# error is (T_h - T_2h) / (2^p - 1)
errors1 = np.abs((comparison2_T1-T_computational2[1:-1]))/(2**p-1)
chiSquared1 = 0
for i in range(len(T_computational2[1:-1])):
    T_observed = comparison2_T1[i]
    x = x_computational2[i+1]
    error = errors1[i]
    T_analytical = analytical_model(x, Pe, Start, End)
    chi_value = (T_observed-T_analytical)**2/error**2
    chiSquared1 += chi_value
reducedChi1 = chiSquared1/len(T_computational2[1:-1])
# The fixed points have no error
errors1 = np.concatenate([[0],errors1,[0]])

errors2 = (np.abs(comparison_T2-T_computational3[1:-1]))/(2**p-1)
chiSquared2 = 0
for i in range(len(T_computational3[1:-1])):
    T_observed = comparison_T2[i]
    x = x_computational3[i+1]
    error = errors2[i]
    T_analytical = analytical_model(x, Pe, Start, End)
    chi_value = (T_observed-T_analytical)**2/error**2
    chiSquared2 += chi_value
reducedChi2 = chiSquared2/len(T_computational3[1:-1])
errors2 = np.concatenate([[0],errors2,[0]])

print(reducedChi1)
print(reducedChi2)
# Formatting preamble
cm = 1/2.54 #cm in inches
mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
# Plotting
fig, axs = plt.subplots(1, 2, figsize=(17*cm, 10*cm), layout="constrained")

axs[0].plot(x_model, y_model, color = "tab:red", linestyle = "-", label = "Exact")
axs[0].errorbar(x_computational2, T_computational1[0::2], yerr=errors1, 
             color = "k", fmt = "x", label = "Calculated")
axs[0].set_title('41 nodes')
axs[0].tick_params(axis='both', direction='in')

axs[1].plot(x_model, y_model, color = "tab:red", linestyle = "-", label = "Exact")
axs[1].errorbar(x_computational3, T_computational2[0::2], yerr=errors2, 
             color = "k", fmt = "x", label = "Calculated")
axs[1].tick_params(axis='both', direction='in')
axs[1].set_title('21 nodes')

handles, labels = axs[1].get_legend_handles_labels()
leg = fig.legend(handles, labels, loc='upper center',
                 bbox_to_anchor=(0.5, 0), ncol = 3)
for ax in axs.flat:
    ax.set(xlabel=r'x/m', ylabel='T/K')
for ax in axs.flat:
    ax.label_outer()
fig.savefig(save_path, bbox_inches = 'tight')
plt.show()
