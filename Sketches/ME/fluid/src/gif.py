import seaborn as sns
import numpy as np
import imageio.v3 as iio
import os
from PIL import Image
rootFolder = os.path.dirname(os.path.dirname(__file__))
dataFolder = os.path.join(rootFolder, "outputs")

specsLocation = os.path.join(dataFolder, "specs.txt")
specsData = np.genfromtxt(specsLocation, delimiter = ",", skip_header = 1)

boxLength = specsData[0]
boxWidth = specsData[1]
holeWidth = specsData[2]
timeStep = specsData[3]
numberOfFiles = specsData[4]

images = []

for i in range(int(numberOfFiles)-5):
    time = i * timeStep
    print(time)
    location  = os.path.join(dataFolder, f"density_t={time:.6f}.txt")
    save_pathpng = os.path.join(dataFolder, f"density_t={time:.1f}.png")
    data = np.genfromtxt(location, delimiter = ",")
    delta = 1
    x = np.arange(0, boxLength, delta)
    y = np.arange(0, boxWidth, delta)
    
    plot = sns.heatmap(data, cmap=sns.color_palette("Blues", as_cmap=True), cbar=False)
    plot.tick_params(left=False, bottom=False)
    plot.set(xticklabels=[])
    plot.set(yticklabels=[])
    fig = plot.get_figure()
    fig.savefig(save_pathpng) 

    images.append(iio.imread(save_pathpng))


pil_frames = [Image.fromarray(f) for f in images]

pil_frames[0].save(os.path.join(dataFolder,'density.gif'), save_all=True, 
                   append_images=pil_frames[1:], duration=300)