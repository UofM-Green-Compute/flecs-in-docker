import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from mpl_toolkits.axisartist.axislines import SubplotZero

def function(x):
    y = 0.5 - 0.25*((x-1)**2)
    return y

def gradient(x, m, c):
    y = m*x + c
    return y

# preamble
cm = 1/2.54 #cm in inches
mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 14})

L = 1.5
h = L/4

root_folder = os.path.dirname(os.path.dirname(__file__))
save_path = os.path.join(root_folder, "outputs", "finite_difference.pdf")

x_list = np.linspace(0,L, 1000)
y_list = function(x_list)

x_points = np.array([h, 2*h, 3*h])
y_points = function(x_points)


m_list = np.array([])
c_list = np.array([])

mUp = (function(2*h)-function(h))/h
cUp = function(h) - mUp*h
m_list = np.append(m_list, mUp)
c_list = np.append(c_list, cUp)

mCentral = (function(3*h)-function(h))/(2*h)
cCentral = function(2*h) - 2*mCentral*h
m_list = np.append(m_list, mCentral)
c_list = np.append(c_list, cCentral)

mDown = (function(3*h)-function(2*h))/h
cDown = function(3*h) - 3*mDown*h
m_list = np.append(m_list, mDown)
c_list = np.append(c_list, cDown)

fig = plt.figure(figsize=(17*cm, 8.5*cm), layout="constrained")
ax = SubplotZero(fig, 1, 1, 1)
fig.add_subplot(ax)
ax.plot(x_list, y_list, color = "tab:red", linewidth=3, zorder=1)
for i in range(len(x_points)):
    ax.scatter(x_points[i], y_points[i], color = "k", s = 100, zorder=3)
    ax.axvline(x_points[i], linestyle = "--", color = "k")
colours = ["tab:blue", "tab:orange", "tab:green"]
gradient_labels = ["Foward", "Central", "Backward"]
for i in range(len(m_list)):
    x_gradient = np.linspace(i*h,(i+2)*h,100)
    y_gradient = gradient(x_gradient, m_list[i], c_list[i])
    ax.plot(x_gradient, y_gradient, color = colours[i], linewidth = 3, zorder = 2, 
            label = gradient_labels[i])
ax.set_ylim(0, 0.7)
ax.set_xlim(0,1.5)
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.tick_params(axis='both', direction='in')
for direction in ["xzero", "yzero"]:
    ax.axis[direction].set_axisline_style("-|>")
    ax.axis[direction].set_visible(True)
    ax.axis[direction].line.set_linewidth(3)
for direction in ["left", "right", "bottom", "top"]:
    ax.axis[direction].set_visible(False)
handles, labels = ax.get_legend_handles_labels()
leg = fig.legend(handles, labels, loc='upper center',
                 bbox_to_anchor=(0.5, 0), ncol = 3)
fig.savefig(save_path, bbox_inches = 'tight')
plt.show()
