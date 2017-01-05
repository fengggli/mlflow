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

    f1 = open(cluster_path, 'r')
    i = 0
    j = 0

    all_classes={}
    prefix = []
    label=[]

    # get all uniq center
    firstline = True

    num_cluster=0

    for line in f1:
        if firstline == True:
            firstline = False
            continue
        words = line.split()

        prefix.append([words[0], words[1], words[2]])
        if words[3] not in all_classes:
            all_classes[words[3]] = num_cluster
            num_cluster += 1
    f1.close()



    f1 = open(cluster_path, 'r')

    count = 0
    firstline = True
    # generate new files
    for line in f1:
        if firstline == True:
            firstline = False
            print("z y x clusterid")
            continue

        words = line.split()
        print(words[0], words[1], words[2], all_classes[words[3]])

    f1.close()




    f1.close()


