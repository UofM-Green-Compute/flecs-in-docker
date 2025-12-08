import os
import numpy as np
import seaborn as sns
import imageio.v3 as iio
from PIL import Image
import matplotlib.pyplot as plt

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

for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    save_pathpdf = os.path.join(dataFolder, f"density_t={time:.1f}.pdf")
    save_pathpng = os.path.join(dataFolder, f"density_t={time:.1f}.png")
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    plot = sns.heatmap(densityData[i], vmin = densityMin, vmax = densityMax, 
                        cmap=sns.color_palette("Blues", as_cmap=True), cbar=False)
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])
    fig = plot.get_figure()
    fig.savefig(save_pathpdf)
    fig.savefig(save_pathpng)
    densityImages.append(iio.imread(save_pathpng))

# Find Horizontal data for each time step and add to array
horizontalData = []
horizontalImages = []
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"horizontal_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    horizontalData.append(data)

horizontalMin, horizontalMax = np.min(horizontalData), np.max(horizontalData)

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

# Find Vertical data for each time step and add to array
verticalData = []
verticalImages = []
for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"vertical_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    verticalData.append(data)
    
verticalMin, verticalMax = np.min(verticalData), np.max(verticalData)

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

vectorImages = []
for i in range(len(horizontalImages)-1):
    time = (i+1)*timeStep
    save_pathpng = os.path.join(dataFolder, f"velocity_field_t={time:.1f}.png")
    save_pathpdf = os.path.join(dataFolder, f"velocity_field_t={time:.1f}.pdf")
    
    u = horizontalData[i+1]
    v = verticalData[i+1]
    x = np.arange(0, boxLength)
    y = np.arange(0, boxWidth)
    X, Y = np.meshgrid(x, y)
    

    plt.figure()
    plt.quiver(X, Y, u, v)
    mask = (u == 0) & (v == 0)
    plt.scatter(X[mask], Y[mask], s=0.01)
    plt.savefig(save_pathpdf)
    plt.savefig(save_pathpng)
    vectorImages.append(iio.imread(save_pathpng))

vectorFrames = [Image.fromarray(f) for f in vectorImages]
vectorFrames[0].save(os.path.join(dataFolder,'velocity_field.gif'), save_all=True, 
                   append_images=vectorFrames[1:], duration=300, loop=0)
