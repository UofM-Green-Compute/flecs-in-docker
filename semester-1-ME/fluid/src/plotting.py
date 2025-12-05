import numpy as np
import seaborn as sns
import os

rootFolder = os.path.dirname(os.path.dirname(__file__))
dataFolder = os.path.join(rootFolder, "outputs")

specsLocation = os.path.join(dataFolder, "specs.txt")
specsData = np.genfromtxt(specsLocation, delimiter = ",", skip_header = 1)

boxLength = specsData[0]
boxWidth = specsData[1]
holeWidth = specsData[2]
timeStep = specsData[3]
numberOfFiles = specsData[4]

# Find density data for each time step and add to array
densityData = []
for i in range(int(numberOfFiles)):
    time = i * timeStep
    location  = os.path.join(dataFolder, f"density_t={time:.6f}.txt")
    save_pathpdf = os.path.join(dataFolder, f"density_t={time:.1f}.pdf")
    data = np.genfromtxt(location, delimiter = ",")
    densityData.append(data)
    print(len(data[0]))
    print(len(data[:,0]))
    delta = 1
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    
    plot = sns.heatmap(data, cmap=sns.color_palette("Blues", as_cmap=True), cbar=False)
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])
    fig = plot.get_figure()
    fig.savefig(save_pathpdf) 

# Find Horizontal data for each time step and add to array
horizontalData = []
for i in range(int(numberOfFiles)):
    time = i * timeStep
    location = specsLocation = os.path.join(dataFolder, f"horizontal_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    horizontalData.append(data)

# Find Vertical data for each time step and add to array
verticalData = []
for i in range(int(numberOfFiles)):
    time = i * timeStep
    location = specsLocation = os.path.join(dataFolder, f"density_t={time:.6f}.txt")
    data = np.genfromtxt(location, delimiter = ",")
    verticalData.append(data)


