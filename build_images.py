import numpy as np
import sys
import matplotlib.pyplot as plt
import multiprocessing
def cell_to_pos(msx, msy, pos):
    return (pos % msx, int(pos / msx))

filename = sys.argv[1]




'''
data  = np.load(filename)
shapex, shapey = data['shape']

print(shapex, shapey)
def processStep(img, i, data=data, shapex=shapex, shapey=shapey):
    image = img
    print('processing step %d'%i)
    grid = np.zeros((shapey, shapex))
    fig, ax = plt.subplots(1, 1, figsize=(20,20))
    for id, type in enumerate(image):
        x, y = cell_to_pos(shapex, shapey, id)
        grid[y][x] = type
    ax.imshow(grid)
    fig.savefig(("{0:0>%d}_water_dummy.jpg" % displ).format(i))
    plt.close(fig)


'''
displ = 4
data = np.load("gids-"+filename)
shapex, shapey = data['shape']
print(data.files)

type = data['type']
print(shapex, shapey)
print(type[0])


def partitionImages(img, i, shapex=shapex, shapey=shapey, type=type):
    image = img
    print('processing %s' % i)

    grid = np.zeros((shapey, shapex))

    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    for id, v in enumerate(image):
        x, y = cell_to_pos(shapex, shapey, id)
        grid[y][x] = v
        print(v)
    print(grid)
    ax.pcolor(grid)
    fig.savefig("{:04d}_partition.jpg".format(int(i.split('_')[-1].split('-')[-1])))
    plt.close(fig)


def processStep(img, i, shapex=shapex, shapey=shapey, type=type):
    image = img
    print('processing %s' % i)

    if(type == 1):
        grid = np.zeros((shapey, shapex))
    else:
        grid = np.ones((shapey, shapex))

    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    for id in image:
        x, y = cell_to_pos(shapex, shapey, id)
        grid[y][x] = type
    ax.imshow(grid)
    fig.savefig("{:04d}_water_dummy.jpg".format(int(i.split('_')[-1].split('-')[-1])))
    plt.close(fig)


cpu_count = multiprocessing.cpu_count()
pool = multiprocessing.Pool(cpu_count)

pool.starmap(partitionImages, [(data[i], i, shapex, shapey, type[0]) for i in [f for f in data.files if 'partition' in f]])
