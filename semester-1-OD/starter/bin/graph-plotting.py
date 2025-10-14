#import matplotlib.pyplot as plt
#import numpy as np

file = open("SHM-Data.txt")

#x = np.linspace(0,100)

pos = []
vel = []
acc = []



for row in file: 
    row = row.split(",")
    pos.append(row[0])
    vel.append(row[1])
    acc.append(row[2]) 

print(pos)

#plt.plot(x,pos)