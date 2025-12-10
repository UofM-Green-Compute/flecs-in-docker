import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import os

def linear_model(x, m, c):
    y = m*x+c
    return y

root_folder = os.path.dirname(os.path.dirname(__file__))
load_path_euler = os.path.join(root_folder, "starter-euler", "outputs", "r2runge.txt")
load_path_runge = os.path.join(root_folder, "starter-runge", "outputs", "r2runge.txt")

euler_r2 = np.genfromtxt(load_path_euler, delimiter=",")
runge_r2 = np.genfromtxt(load_path_runge, delimiter=",")

logTimeEuler = np.log(euler_r2[:,0])
logR2Euler = np.log(1-euler_r2[:,1])

logTimeRunge = np.log(runge_r2[:,0])
logR2Runge = np.log(1-runge_r2[:,1])

poptEuler, pcovEuler = curve_fit(linear_model, logTimeEuler, logR2Euler)
poptRunge, pcovRunge = curve_fit(linear_model, logTimeRunge, logR2Runge)

x_runge = np.linspace(np.min(logTimeRunge), np.max(logTimeRunge), 1000)
x_euler = np.linspace(np.min(logTimeEuler), np.max(logTimeEuler), 1000)
euler_model = linear_model(x_euler, poptEuler[0], poptEuler[1])
runge_model = linear_model(x_runge, poptRunge[0], poptRunge[1])

print(poptEuler[0], np.e**poptEuler[1])
print(poptRunge[0], np.e**poptRunge[1])

plt.figure()
plt.scatter(logTimeEuler, logR2Euler, marker = "x", label = "Euler Method")
plt.scatter(logTimeRunge, logR2Runge, marker = "x", label = "Runge-Kutta Method")
plt.plot(x_euler, euler_model, label = "euler fit")
plt.plot(x_runge, runge_model, label = "runge fit")
plt.xlabel(r"ln($\Delta$t)")
plt.ylabel(r"ln(1-$r^2$)")
plt.legend()
plt.show()