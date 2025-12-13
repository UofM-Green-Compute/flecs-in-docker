import os
import numpy as np
import seaborn as sns
import imageio.v3 as iio
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl


# Formatting preamble
cm = 1/2.54 #cm in inches
mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 12})

rootFolder = os.path.dirname(os.path.dirname(__file__))
dataFolder = os.path.join(rootFolder, "outputs")

specsLocation = os.path.join(dataFolder, "specs.txt")
specsData = np.genfromtxt(specsLocation, delimiter = ",", skip_header = 1)

boxLength = specsData[0]
boxWidth = specsData[1]
holeWidth = specsData[2]
timeStep = specsData[3]
numberOfFiles = specsData[4]

delta = 1

# Find density data for each time step and add to array
densityData = []
densityImages = []
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"density_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    densityData.append(data)

densityMin, densityMax = np.min(densityData), np.max(densityData)

# Find Horizontal data for each time step and add to array
horizontalData = []
horizontalImages = []
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"horizontal_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    horizontalData.append(data)

horizontalMin, horizontalMax = np.min(horizontalData), np.max(horizontalData)

# Find Vertical data for each time step and add to array
verticalData = []
verticalImages = []
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"vertical_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    verticalData.append(data)
    
verticalMin, verticalMax = np.min(verticalData), np.max(verticalData)

timePlots = 4
final_save_path = os.path.join(dataFolder, f"report_fluids.pdf")
densityPlot = sns.heatmap(densityData[i], vmin=densityMin, vmax=densityMax, cmap='gray_r')
horizontalData = np.array([horizontalData])
verticalData = np.array([verticalData])
speedData = np.sqrt(horizontalData**2+verticalData**2)
speedMin = np.min(speedData)
speedMax = np.max(speedData)
titles = ["Density", "Streamfunction", "Fluid Speed"]
fig, axs = plt.subplots(4, 3, sharex='col', figsize=(17*cm, 12*cm), layout="constrained")
for i in range(timePlots):
    for j in range(3):
        axs[i,j].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
        if (i == 0):
            axs[0, j].set_title(titles[j])
    u = horizontalData[0,i+1]
    v = verticalData[0,i+1]
    x = np.arange(0, boxLength)
    y = np.arange(0, boxWidth)
    X, Y = np.meshgrid(x, y)
    densityPlot = sns.heatmap(densityData[i+1], vmin=densityMin, vmax=densityMax, cmap = 'gray_r', 
                              ax = axs[i, 0], cbar = False)
    for _, spine in densityPlot.spines.items():
        spine.set_visible(True)
    axs[i,1].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    speedPlot = sns.heatmap(speedData[0,i+1], vmin=speedMin, vmax = speedMax, cmap = 'gray_r', 
                            ax = axs[i, 2], cbar = False)
    for _, spine in speedPlot.spines.items():
        spine.set_visible(True)
fig.savefig(final_save_path, bbox_inches = 'tight')
plt.show()

'''
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    save_pathpdf = os.path.join(dataFolder, f"density_t={time:.1f}.pdf")
    save_pathpng = os.path.join(dataFolder, f"density_t={time:.1f}.png")
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    rows, cols = densityData[i].shape

    # Make figure size proportional to data shape
    plt.figure(figsize=(17*cm, 8.5*cm))  # or multiply both by any constant

    plot = sns.heatmap(densityData[i], vmin=densityMin, vmax=densityMax, cmap='gray_r')
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])

    plt.tight_layout(pad=0)
    fig = plot.get_figure()
    fig.savefig(save_pathpdf, bbox_inches='tight', pad_inches=0)
    fig.savefig(save_pathpng, bbox_inches='tight', pad_inches=0)
    densityImages.append(iio.imread(save_pathpng))

for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    save_pathpdf = os.path.join(dataFolder, f"horizontal_t={time:.1f}.pdf")
    save_pathpng = os.path.join(dataFolder, f"horizontal_t={time:.1f}.png")
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    plot = sns.heatmap(horizontalData[i], vmin = horizontalMin, vmax = horizontalMax, 
                        cmap=sns.color_palette("Blues", as_cmap=True), cbar=False)
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])
    fig = plot.get_figure()
    fig.savefig(save_pathpdf) 
    fig.savefig(save_pathpng)
    horizontalImages.append(iio.imread(save_pathpng))

for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    save_pathpdf = os.path.join(dataFolder, f"vertical_t={time:.1f}.pdf")
    save_pathpng = os.path.join(dataFolder, f"vertical_t={time:.1f}.png")
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    plot = sns.heatmap(verticalData[i], vmin = verticalMin, vmax = verticalMax, 
                       cmap=sns.color_palette("Blues", as_cmap=True), cbar=False)
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])
    fig = plot.get_figure()
    fig.savefig(save_pathpdf) 
    fig.savefig(save_pathpng) 
    verticalImages.append(iio.imread(save_pathpng))

densityFrames = [Image.fromarray(f) for f in densityImages]
densityFrames[0].save(os.path.join(dataFolder,'density.gif'), save_all=True, 
                   append_images=densityFrames[1:], duration=300, loop=0)

horizontalFrames = [Image.fromarray(f) for f in horizontalImages]
horizontalFrames[0].save(os.path.join(dataFolder,'horizontal_velocity.gif'), save_all=True, 
                   append_images=densityFrames[1:], duration=300, loop=0)

verticalFrames = [Image.fromarray(f) for f in verticalImages]
verticalFrames[0].save(os.path.join(dataFolder,'vertical_velocity.gif'), save_all=True, 
                   append_images=verticalFrames[1:], duration=300, loop=0)

quiverImages = []
streamImages = []
for i in range(len(horizontalImages)-1):
    time = (i+1)*timeStep
    save_pathpng1 = os.path.join(dataFolder, f"velocity_field_t={time:.1f}.png")
    save_pathpdf1 = os.path.join(dataFolder, f"velocity_field_t={time:.1f}.pdf")
    save_pathpng2 = os.path.join(dataFolder, f"streamfunction_t={time:.1f}.png")
    save_pathpdf2 = os.path.join(dataFolder, f"streamfunction_t={time:.1f}.pdf")
    
    u = horizontalData[i+1]
    v = verticalData[i+1]
    x = np.arange(0, boxLength)
    y = np.arange(0, boxWidth)
    X, Y = np.meshgrid(x, y)

    plt.figure(figsize=(17*cm, 8.5*cm)) 
    plt.tight_layout(pad=0)
    plt.quiver(X[::3,::3], Y[::3,::3], u[::3,::3], v[::3,::3])
    plt.savefig(save_pathpdf1)
    plt.savefig(save_pathpng1)
    quiverImages.append(iio.imread(save_pathpng1))

    plt.figure(figsize=(17*cm, 8.5*cm)) 
    plt.tight_layout(pad=0)
    plt.streamplot(X, Y, u, v, color = 'k', density = 3)
    plt.savefig(save_pathpdf2)
    plt.savefig(save_pathpng2)
    streamImages.append(iio.imread(save_pathpng2))

quiverFrames = [Image.fromarray(f) for f in quiverImages]
quiverFrames[0].save(os.path.join(dataFolder,'velocity_field.gif'), save_all=True, 
                   append_images=quiverFrames[1:], duration=300, loop=0)

streamFrames = [Image.fromarray(f) for f in streamImages]
streamFrames[0].save(os.path.join(dataFolder,'streamfuntion.gif'), save_all=True, 
                   append_images=streamFrames[1:], duration=300, loop=0)
'''