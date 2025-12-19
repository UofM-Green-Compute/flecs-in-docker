import os
import numpy as np
import seaborn as sns
import imageio.v3 as iio
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import ScaledTranslation

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
nodeSpacing = 0.01

delta = 1

# Find density data for each time step and add to array
densityData = []
densityImages = []
for i in range(int(numberOfFiles)-6):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"density_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    densityData.append(data)

densityMin, densityMax = np.min(densityData), np.max(densityData)

# Find Horizontal data for each time step and add to array
horizontalData = []
horizontalImages = []
for i in range(int(numberOfFiles)-6):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"horizontal_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    horizontalData.append(data)

horizontalMin, horizontalMax = np.min(horizontalData), np.max(horizontalData)

# Find Vertical data for each time step and add to array
verticalData = []
verticalImages = []
for i in range(int(numberOfFiles)-6):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"vertical_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    verticalData.append(data)
    
verticalMin, verticalMax = np.min(verticalData), np.max(verticalData)

timePlots = 4
final_save_path1 = os.path.join(dataFolder, f"fluids_density.pdf")
final_save_pathpng = os.path.join(dataFolder, f"fluids_density.png")
horizontalData = np.array([horizontalData])
verticalData = np.array([verticalData])
speedData = np.sqrt(horizontalData**2+verticalData**2)
speedMin = np.min(speedData)
speedMax = np.max(speedData)
fig, axs = plt.subplot_mosaic([['A)', 'B)'], ['C)', 'D)']], sharey=True, sharex=True, figsize=(17*cm, 10*cm), layout="constrained")
for i in range(4):
    x = np.arange(0, boxLength*100+1, dtype=int)
    y = np.arange(0, boxWidth*100+1, dtype = int)
    z = densityData[i+1]
    if i == 0:
        c1=axs['A)'].pcolormesh(x, y, z, cmap='Purples', vmin=densityMin, vmax=densityMax)
    if i == 1:
        c2=axs['B)'].pcolormesh(x, y, z, cmap='Purples', vmin=densityMin, vmax=densityMax)
    if i == 2:
        c3=axs['C)'].pcolormesh(x, y, z, cmap='Purples', vmin=densityMin, vmax=densityMax)
    if i == 3:
        c4=axs['D)'].pcolormesh(x, y, z, cmap='Purples', vmin=densityMin, vmax=densityMax)

axs['C)'].set_xlabel(r'x/cm')
axs['D)'].set_xlabel(r'x/cm')

axs['A)'].set_ylabel(r'y/cm')
axs['C)'].set_ylabel(r'y/cm')

axs['A)'].tick_params(axis='both', direction='in')
axs['B)'].tick_params(axis='both', direction='in')
axs['C)'].tick_params(axis='both', direction='in')
axs['D)'].tick_params(axis='both', direction='in')

axs['A)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['B)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['C)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['D)'].axis([x.min(), x.max(), y.min(), y.max()])

for label, ax in axs.items():
    ax.text( 0.0, 1.0, label, 
            transform=(ax.transAxes + ScaledTranslation(+5/72, -5/72, fig.dpi_scale_trans)),
            fontsize='medium', va='top', ha='left', color='white'
    )
    
fig.colorbar(c1, ax=list(axs.values()))
fig.savefig(final_save_path1, bbox_inches = 'tight')
fig.savefig(final_save_pathpng, dpi = 600, bbox_inches = 'tight')
plt.show()

final_save_path2 = os.path.join(dataFolder, f"speed.pdf")
fig, axs = plt.subplot_mosaic([['A)', 'B)'], ['C)', 'D)']], sharey=True, sharex=True, figsize=(17*cm, 10*cm), layout="constrained")
for i in range(4):
    x = np.arange(0, boxLength*100+1, dtype=int)
    y = np.arange(0, boxWidth*100+1, dtype = int)
    z = speedData[0,i+1]
    if i == 0:
        c1=axs['A)'].pcolormesh(x, y, z, cmap='Reds', vmin=speedMin, vmax=speedMax)
    if i == 1:
        c2=axs['B)'].pcolormesh(x, y, z, cmap='Reds', vmin=speedMin, vmax=speedMax)
    if i == 2:
        c3=axs['C)'].pcolormesh(x, y, z, cmap='Reds', vmin=speedMin, vmax=speedMax)
    if i == 3:
        c4=axs['D)'].pcolormesh(x, y, z, cmap='Reds', vmin=speedMin, vmax=speedMax)

axs['C)'].set_xlabel(r'x/cm')
axs['D)'].set_xlabel(r'x/cm')

axs['A)'].set_ylabel(r'y/cm')
axs['C)'].set_ylabel(r'y/cm')

axs['A)'].tick_params(axis='both', direction='in')
axs['B)'].tick_params(axis='both', direction='in')
axs['C)'].tick_params(axis='both', direction='in')
axs['D)'].tick_params(axis='both', direction='in')

axs['A)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['B)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['C)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['D)'].axis([x.min(), x.max(), y.min(), y.max()])

for label, ax in axs.items():
    ax.text( 0.0, 1.0, label, 
            transform=(ax.transAxes + ScaledTranslation(+5/72, -5/72, fig.dpi_scale_trans)),
            fontsize='medium', va='top', ha='left'
    ) 
fig.colorbar(c1, ax=list(axs.values()))
fig.savefig(final_save_path2, dpi = 600, bbox_inches = 'tight')
plt.show()

final_save_path3 = os.path.join(dataFolder, f"streamfunctions.pdf")
fig, axs = plt.subplot_mosaic([['A)', 'B)'], ['C)', 'D)']], sharey=True, sharex=True, figsize=(17*cm, 10*cm), layout="constrained")
for i in range(4):
    u = horizontalData[0,i+1]
    v = verticalData[0,i+1]
    x = np.arange(0, boxLength*100)
    y = np.arange(0, boxWidth*100)
    X, Y = np.meshgrid(x, y)
    if i == 0:
        axs['A)'].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 1:
        axs['B)'].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 2:
        axs['C)'].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    if i == 3:
        axs['D)'].streamplot(X, Y, u, v, color = 'k', density = 2, linewidth=0.1)
    
axs['C)'].set_xlabel(r'x/cm')
axs['D)'].set_xlabel(r'x/cm')

axs['A)'].set_ylabel(r'y/cm')
axs['C)'].set_ylabel(r'y/cm')

axs['A)'].tick_params(axis='both', direction='in')
axs['B)'].tick_params(axis='both', direction='in')
axs['C)'].tick_params(axis='both', direction='in')
axs['D)'].tick_params(axis='both', direction='in')

axs['A)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['B)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['C)'].axis([x.min(), x.max(), y.min(), y.max()])
axs['D)'].axis([x.min(), x.max(), y.min(), y.max()])

for label, ax in axs.items():
    ax.text( 0.0, 1.0, label, 
            transform=(ax.transAxes + ScaledTranslation(+5/72, -5/72, fig.dpi_scale_trans)),
            fontsize='medium', va='top', ha='left'
    ) 
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