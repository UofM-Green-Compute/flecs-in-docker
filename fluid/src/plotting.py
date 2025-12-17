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
    print(f'{time:.6f}')
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
final_save_path1 = os.path.join(dataFolder, f"fluids_density.pdf")
horizontalData = np.array([horizontalData])
verticalData = np.array([verticalData])
speedData = np.sqrt(horizontalData**2+verticalData**2)
speedMin = np.min(speedData)
speedMax = np.max(speedData)

fig, axs = plt.subplots(2, 2, sharey='row', figsize=(17*cm, 10*cm), layout="constrained")
for i in range(4):
    print(i)
    print(densityData[i+1])
    if i == 0:
        densityPlot1 = sns.heatmap(densityData[i+1], vmin=densityMin, vmax=densityMax, cmap = sns.color_palette("flare", as_cmap=True), 
                              ax = axs[0, 0], cbar = False)
        for _, spine in densityPlot1.spines.items():
            spine.set_visible(True)
    if i == 1:
        densityPlot2 = sns.heatmap(densityData[i+1], vmin=densityMin, vmax=densityMax, cmap = sns.color_palette("flare", as_cmap=True), 
                              ax = axs[0, 1], cbar = False)
        for _, spine in densityPlot2.spines.items():
            spine.set_visible(True)
    if i == 2:
        densityPlot3 = sns.heatmap(densityData[i+1], vmin=densityMin, vmax=densityMax, cmap = sns.color_palette("flare", as_cmap=True), 
                              ax = axs[1, 0], cbar = False)
        for _, spine in densityPlot3.spines.items():
            spine.set_visible(True)
    if i == 3:
        densityPlot4 = sns.heatmap(densityData[i+1], vmin=densityMin, vmax=densityMax, cmap = sns.color_palette("flare", as_cmap=True), 
                              ax = axs[1, 1], cbar = False)
        for _, spine in densityPlot4.spines.items():
            spine.set_visible(True)
    axs[0,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[0,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    
fig.suptitle('Density', fontsize=20)
mappable1 = densityPlot1.collections[0]
fig.colorbar(mappable1, ax=axs, orientation='vertical', location="right", fraction=.1)
fig.savefig(final_save_path1, dpi = 600, bbox_inches = 'tight')
plt.show()

final_save_path2 = os.path.join(dataFolder, f"speed.pdf")
fig, axs = plt.subplots(2, 2, sharex='col', figsize=(17*cm, 10*cm), layout="constrained")
for i in range(4):
    u = horizontalData[0,i+1]
    v = verticalData[0,i+1]
    x = np.arange(0, boxLength)
    y = np.arange(0, boxWidth)
    X, Y = np.meshgrid(x, y)
    if i == 0:
        speedPlot1 = sns.heatmap(speedData[0,i+1], vmin=speedMin, vmax = speedMax, 
                                cmap = sns.color_palette("YlOrBr", as_cmap=True), 
                                ax = axs[0, 0], cbar = False)
        for _, spine in speedPlot1.spines.items():
            spine.set_visible(True)
    if i == 1:
        speedPlot2 = sns.heatmap(speedData[0,i+1], vmin=speedMin, vmax = speedMax, 
                                cmap = sns.color_palette("YlOrBr", as_cmap=True), ax = axs[0, 1], 
                                cbar = False)
        for _, spine in speedPlot2.spines.items():
            spine.set_visible(True)
    if i == 2:
        speedPlot3 = sns.heatmap(speedData[0,i+1], vmin=speedMin, vmax = speedMax, 
                                cmap = sns.color_palette("YlOrBr", as_cmap=True), ax = axs[1, 0],
                                cbar = False)
        for _, spine in speedPlot3.spines.items():
            spine.set_visible(True)
    if i == 3:
        speedPlot4 = sns.heatmap(speedData[0,i+1], vmin=speedMin, vmax = speedMax, 
                                cmap = sns.color_palette("YlOrBr", as_cmap=True), ax = axs[1, 1],
                                cbar = False)
        for _, spine in speedPlot4.spines.items():
            spine.set_visible(True)
    axs[0,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[0,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    

axs[1,0].set_xticks([0, 50, 100, 150, 200])
axs[1,1].set_xticks([0, 50, 100, 150, 200])

axs[1,0].set_yticks([0, 50, 100])
axs[0,0].set_yticks([0, 50, 100])

axs[1,0].set_ylabel(r'y/cm')
axs[0,0].set_ylabel(r'y/cm')
 
mappable2 = speedPlot1.collections[0]
fig.suptitle('Speed', fontsize=20)
fig.colorbar(mappable1, ax=axs, orientation='vertical', location="right", fraction=.1)
fig.savefig(final_save_path2, dpi = 600, bbox_inches = 'tight')
plt.show()

final_save_path3 = os.path.join(dataFolder, f"streamfunctions.pdf")
fig, axs = plt.subplots(2, 2, sharex='col', figsize=(17*cm, 12*cm), layout="constrained")
for i in range(4):
    u = horizontalData[0,i+1]
    v = verticalData[0,i+1]
    x = np.arange(0, boxLength)
    y = np.arange(0, boxWidth)
    X, Y = np.meshgrid(x, y)
    if i == 0:
        axs[0,0].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 1:
        axs[0,1].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 2:
        axs[1,0].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 3:
        axs[1,1].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    axs[0,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[0,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,0].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )
    axs[1,1].tick_params(
            bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )

fig.suptitle('Streamfunction', fontsize=20)
fig.savefig(final_save_path3, bbox_inches = 'tight')
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