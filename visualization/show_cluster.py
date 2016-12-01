from matplotlib.colors import ListedColormap, NoNorm
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys

if(len(sys.argv)<2):
    print("python show_cluster.py clusterid_text_file")
    quit()
else:

    cluster_path=sys.argv[1]

    # show clustering result
    sideLength = 20
    matrix = np.zeros((sideLength, sideLength), dtype=int)
    f1 = open(cluster_path, 'r')
    i = 0
    j = 0

    all_classes=[]
    for line in f1:
        tmp = int(line)
        if tmp not in all_classes:
            all_classes.append(tmp)
    f1.close();


    f1 = open(cluster_path, 'r')
    print(all_classes)
    for line in f1:
        tmp = int(line)
        if tmp == all_classes[0]:
            matrix[i,j] = 1
        elif tmp == all_classes[1]:
            matrix[i,j] = 2
        else:
            matrix[i,j] = 3

        j = j + 1
        if j == sideLength:
            i = i+1
            j = 0
    f1.close();

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

    img = plt.imread('visualization/nogrid.PNG')
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
