import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib
from MyDelaunay import MyDelaunay
import time
import os.path

#####Read in data -- already downsampled #####
input_file_name = 'BoulderFalls_1000.csv' ##File name with .csv extension
with open(input_file_name) as f:
    next(csv.reader(f)) ##skip header
    data = list(csv.reader(f, quoting = csv.QUOTE_NONNUMERIC))
  
##Data 
a = np.asarray(data)

##Triangle object
tri = MyDelaunay(np.vstack([a[:,0],a[:,1]]).T, pts = a)


##Source of sound
source = np.asarray([465100,  4.42818E6, 2130])


##Intensity of source
source_intensity = 100 # decibels of input source
d = tri.distance(source) # distance from source to each simplex
inv_mag = np.divide(1,np.power(d,2))

#The step_factor can be set to speed up computation time
##Source to destination blocking vector
source_blocked = tri.check_blocked(source, step_factor = 20)

#Element-wise initial intensity due to source
#First term is sound propogation inverse square law
I0 = np.multiply(source_intensity - 20 * np.log10(d), source_blocked)

##Check if the element to element blocking matrix has already been 
##created. If so, simply read it in. Otherwise build it
if os.path.exists('block_mat_' + input_file_name):
    with open('block_mat_' + input_file_name) as f:
        data_bm = list(csv.reader(f, quoting = csv.QUOTE_NONNUMERIC))
    block_mat = np.asarray(data_bm)
    
else:    
    block_mat = tri.element_blocking_matrix()
    np.savetxt('block_mat_' + input_file_name ,block_mat, delimiter=',', fmt = '%.0f')


##Include first order reflections?
ref_bool = True
if ref_bool:
    ref = tri.reflection_intensity(source, block_mat, I_incident = I0, dispersion_angle = 0.4)
    I = I0 + ref  #reflection + I0
else:
    I = I0


# %% Plotting
ax = plt.axes(projection='3d')

#Decibel scales
loud = 85
quiet = 40
norm = matplotlib.colors.LogNorm(quiet, loud)
colors = matplotlib.cm.jet(norm(I))

triang = mtri.Triangulation(x=a[:, 0], y=a[:, 1], triangles=tri.simplices)
surf = ax.plot_trisurf(triang, a[:,2], ax.view_init(30, 30),
                edgecolor='black', 
                linewidth = 0.2)
surf.set_fc(colors)

plt.rcParams['figure.dpi'] = 300
plt.show()

 