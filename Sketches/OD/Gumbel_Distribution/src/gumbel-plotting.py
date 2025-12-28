import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Open data file
file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_Dust.txt")
specs_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_specifications.txt"); 
density_file = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/density_field.txt")
density_gambel_20 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=20.txt")
density_gambel_50 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=50.txt")
density_gambel_75 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=75.txt")
density_gambel_100 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=100.txt")
density_gambel_200 = open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/Gambel_Density_ParticleNo=200.txt")

# *** READ IN DATA ***

# Read Specifications file - info about the constants and parameters in the code 
for i,row in enumerate(specs_file): 
    if (i == 0):
        row = row.split("|")
        DT = float(row[0].split(":")[1])
        no_steps = int(row[1].split(":")[1])
        no_particles = int(row[2].split(":")[1])
    if (i == 1):
        row = row.split("|") 
        Temp = int(row[0].split(":")[1])
        Mass = float(row[1].split(":")[1])
    if (i == 2):
        row = row.split("|") 
        GX = int(row[0].split(":")[1])
        GX2 = int(row[1].split(":")[1])
        GY = int(row[2].split(":")[1])
        slit_width = int(row[3].split(":")[1])
    if (i == 3):
        row = row.split("|") 
        no_gambel_points = int(row[0].split(":")[1])
        no_gambel_averages = int(row[1].split(":")[1])

gambel_matrix_20 = np.zeros((no_gambel_points,2))
gambel_matrix_50 = np.zeros((no_gambel_points,2))
gambel_matrix_75 = np.zeros((no_gambel_points,2))
gambel_matrix_100 = np.zeros((no_gambel_points,2))
gambel_matrix_200 = np.zeros((no_gambel_points,2))

# Read in data on how estimation error varies with h
count20 = 1
factor20 = 0 
for i,row in enumerate(density_gambel_20): 
    row = row.split(",")
    for j,rho in enumerate(row):
        gambel_matrix_20[i-no_gambel_points*factor20][j] += float(rho) / no_gambel_averages
    count20 += 1
    if(count20 == no_gambel_points):
        count20 = 0
        factor20 += 1


count50 = 1
factor50 = 0
for i,row in enumerate(density_gambel_50): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_50[i-no_gambel_points*factor50][j] += float(rho) / no_gambel_averages
    count50 += 1
    if(count50 == no_gambel_points):
        count50 = 0
        factor50 += 1

count75 = 1
factor75 = 0
for i,row in enumerate(density_gambel_75): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_75[i-no_gambel_points*factor75][j] += float(rho) / no_gambel_averages
    count75 += 1
    if(count75 == no_gambel_points):
        count75 = 0
        factor75 += 1

count100 = 1
factor100 = 0
for i,row in enumerate(density_gambel_100): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_100[i-no_gambel_points*factor100][j] += float(rho) / no_gambel_averages
    count100 += 1
    if(count100 == no_gambel_points):
        count100 = 0
        factor100 += 1

count200 = 1
factor200 = 0
for i,row in enumerate(density_gambel_200): 
    row = row.split(",")
    for j,rho in enumerate(row):
            gambel_matrix_200[i-no_gambel_points*factor200][j] += float(rho) / no_gambel_averages
    count200 += 1
    if(count200 == no_gambel_points):
        count200 = 0
        factor200 += 1

# *** FORMAT GRAPHS ***

# Create figures and axes 
cm = 1/2.54 #cm in inches
fig4, ax_gambel = plt.subplots(figsize=(17*cm,10*cm))

# Format graphs
mpl.rcParams.update({'font.size': 12})

# Major grid:
ax_gambel.grid(True, which='major', linestyle='-', linewidth=0.75, alpha=0.25)
# Minor ticks and grid:
ax_gambel.minorticks_on()
ax_gambel.grid(True, which='minor', linestyle='-', linewidth=0.25, alpha=0.15)
ax_gambel.set_axisbelow(True) # <-- Ensure grid is below data
ax_gambel.tick_params(axis='both',direction='in')

ax_gambel.set_xlabel("h")
ax_gambel.set_ylabel("Average Error")

# Scatter plot of Gambel errors 
ax_gambel.scatter(gambel_matrix_20[:,0], gambel_matrix_20[:,1], s=15, color='k', marker='.',label = '20 Samples')
ax_gambel.scatter(gambel_matrix_50[:,0], gambel_matrix_50[:,1], s=15, color='k', marker='^',label = '50 Samples')
ax_gambel.scatter(gambel_matrix_75[:,0], gambel_matrix_75[:,1], s=15, color='k', marker='>',label = '75 Samples')
ax_gambel.scatter(gambel_matrix_100[:,0], gambel_matrix_100[:,1], s=15, color='k', marker=',',label = '100 Samples')
ax_gambel.scatter(gambel_matrix_200[:,0], gambel_matrix_200[:,1], color='k', marker='o',label = '200 Samples')

plt.show()