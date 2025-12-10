import matplotlib.pyplot as plt
import numpy as np
import os

root_folder = os.path.dirname(os.path.dirname(__file__))
load_path_euler = os.path.join(root_folder, "starter-euler", "outputs", "r2runge.txt")
load_path_runge = os.path.join(root_folder, "starter-runge", "outputs", "r2runge.txt")

euler_r2 = np.genfromtxt(load_path_euler, delimiter=",")
runge_r2 = np.genfromtxt(load_path_runge, delimiter=",")

plt.figure()
plt.scatter(np.log(euler_r2[:,0]), np.log(1-euler_r2[:,1]), marker = "x", label = "Euler Method")
plt.scatter(np.log(runge_r2[:,0]), np.log(1-runge_r2[:,1]), marker = "x", label = "Runge-Kutta Method")
plt.xlabel(r"ln($\Delta$t)")
plt.ylabel(r"ln(1-$r^2$)")
plt.legend()
plt.show()