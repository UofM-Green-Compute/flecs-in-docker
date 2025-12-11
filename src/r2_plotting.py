import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
import numpy as np
import os

def linear_model(x, m, c):
    y = m*x+c
    return y

root_folder = os.path.dirname(os.path.dirname(__file__))
load_path_eulerR2 = os.path.join(root_folder, "starter-euler", "outputs", "r2runge.txt")
load_path_rungeR2 = os.path.join(root_folder, "starter-runge", "outputs", "r2runge.txt")
load_path_eulerEnergy = os.path.join(root_folder, "starter-euler", "outputs", "energy.txt")
load_path_rungeEnergy = os.path.join(root_folder, "starter-runge", "outputs", "energy.txt")
save_path = os.path.join(root_folder, "outputs", "oscillator_method_analysis.pdf")

eulerDataR2 = np.genfromtxt(load_path_eulerR2, delimiter=",")
rungeDataR2 = np.genfromtxt(load_path_rungeR2, delimiter=",")

eulerDataEnergy = np.genfromtxt(load_path_eulerEnergy, delimiter=",")
rungeDataEnergy = np.genfromtxt(load_path_rungeEnergy, delimiter=",")

TimeEuler = eulerDataR2[:,0]
logTimeEuler = np.log(TimeEuler)
r2Euler = eulerDataR2[:,1]
logR2Euler = np.log(r2Euler)
logTransformedEuler = np.log(1-r2Euler)
EnergyEuler = eulerDataEnergy[:,1]
logEnergyEuler = np.log(EnergyEuler)

TimeRunge = rungeDataR2[:,0]
logTimeRunge = np.log(TimeRunge)
r2Runge = rungeDataR2[:,1]
logR2Runge = np.log(r2Runge)
logTransformedRunge = np.log(1-r2Runge)
EnergyRunge = rungeDataEnergy[:,1]
logEnergyRunge = np.log(EnergyRunge)

poptEulerR2, pcovEulerR2 = curve_fit(linear_model, logTimeEuler, logTransformedEuler)
poptRungeR2, pcovRungeR2 = curve_fit(linear_model, logTimeRunge, logTransformedRunge)
poptEulerEnergy, pcovEulerEnergy = curve_fit(linear_model, logTimeEuler[2:], logEnergyEuler[2:])
poptRungeEnergy, pcovRungeEnergy = curve_fit(linear_model, logTimeRunge[1:], logEnergyRunge[1:])

logtR2Runge = np.linspace(np.min(logTimeRunge), np.max(logTimeRunge), 1000)
logtR2Euler = np.linspace(np.min(logTimeEuler), np.max(logTimeEuler), 1000)
eulerR2model = linear_model(logtR2Euler, poptEulerR2[0], poptEulerR2[1])
rungeR2model = linear_model(logtR2Runge, poptRungeR2[0], poptRungeR2[1])
logtEnergyRunge = np.linspace(np.min(logTimeRunge[1:]), np.max(logTimeRunge[1:]), 1000)
logtEnergyEuler = np.linspace(np.min(logTimeEuler[2:]), np.max(logTimeEuler[2:]), 1000)
eulerEnergymodel = linear_model(logtEnergyEuler, poptEulerEnergy[0], poptEulerEnergy[1])
rungeEnergymodel = linear_model(logtEnergyRunge, poptRungeEnergy[0], poptRungeEnergy[1])

# prints all the fitting parameters
print(poptEulerR2[0], np.e**poptEulerR2[1])
print(poptRungeR2[0], np.e**poptRungeR2[1])
print(-poptEulerEnergy[0], np.e**poptEulerEnergy[1])
print(-poptRungeEnergy[0], np.e**poptRungeEnergy[1])

# Formatting preamble
cm = 1/2.54 #cm in inches
mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})
# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17*cm, 10*cm), layout="constrained")
ax1.scatter(logTimeEuler, logTransformedEuler, marker = "x", label = "Euler Method")
ax1.scatter(logTimeRunge, logTransformedRunge, marker = "x", label = "Runge-Kutta Method")
ax1.plot(logtR2Euler, eulerR2model, label = "euler fit")
ax1.plot(logtR2Runge, rungeR2model, label = "runge fit")
ax1.set_xlabel(r"ln($\Delta$t/s)")
ax1.set_ylabel(r"ln(1-$r^2$)")
ax1.tick_params(axis='both', direction='in')


ax2.scatter(logTimeEuler[2:], logEnergyEuler[2:], marker = "x", color = "tab:blue", s=20,
            label = "Euler Method")
ax2.scatter(logTimeRunge[1:], logEnergyRunge[1:], marker = "x", color = "tab:orange", s = 20,
            label = "Runge-Kutta Method")
ax2.scatter(logTimeEuler[:2], logEnergyEuler[:2], marker = "s", color = "tab:blue", s = 20,
            label = "Euler Anomalies")
ax2.scatter(logTimeRunge[:1], logEnergyRunge[:1], marker = "s", color = "tab:orange", s = 20, 
            label = "Runge-Kutta Anomalies")
ax2.plot(logtEnergyEuler, eulerEnergymodel, label = "euler fit")
ax2.plot(logtEnergyRunge, rungeEnergymodel, label = "runge fit")
ax2.set_xlabel(r"ln($\Delta$t/s)")
ax2.set_ylabel(r"ln(Energy/J)")
ax2.tick_params(axis='both', direction='in')

handles, labels = ax2.get_legend_handles_labels()
leg = fig.legend(handles, labels, loc='upper center',
                 bbox_to_anchor=(0.5, 0), ncol = 3)

fig.savefig(save_path, bbox_inches = 'tight')
plt.show()
