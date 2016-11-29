from matplotlib.colors import ListedColormap, NoNorm
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# show clustering result
sideLength = 20
matrix = np.zeros((sideLength, sideLength), dtype=int)
f1 = open('data/clusterid.txt', 'r')
i = 0
j = 0
for line in f1:
    tmp = int(line)
    if tmp == 194:
        matrix[i,j] = 1
    elif tmp == 309:
        matrix[i,j] = 3
    else:
        matrix[i,j] = 2

    j = j + 1
    if j == sideLength:
        i = i+1
        j = 0
cmap = ListedColormap(['#E0E0E0', '#FF8C00', '#8c00FF', '#00FF8C'])

gs = gridspec.GridSpec(1,2,width_ratios=[1,2],
                       height_ratios=[2,2])
plt.subplot(gs[0,0])

plt.pcolor(matrix,cmap=cmap,norm=NoNorm(),)
plt.grid(True, which='both', ls='-')
plt.title('clustering results')
ax = plt.gca()
ax.invert_yaxis()
ax.set_xticks(range(0,20))
ax.set_yticks(range(0,20))

plt.xlim(0,20)
#plt.show()

# also show real data in paraview
plt.subplot(gs[0,1])

img = plt.imread('data/nogrid.PNG')
plt.title('real data')
#modyfy the image
#img[:,::dy,:] = grid_color
#img[::dx,:,:] = grid_color
plt.imshow(img, extent=[0,20,20,0])
plt.grid(True, which='both',ls='-')
ax=plt.gca()
ax.set_xticks(range(0,20))
ax.set_yticks(range(0,20))
#img.grid(True, color='r', linestyle='--', linewidth=2)
plt.xlim(0,20)
plt.tight_layout()
plt.show()
