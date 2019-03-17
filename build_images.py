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
    fig.savefig(("{0:0>%d}_water_dummy.pdf" % displ).format(i.split('_')[-1]))
    plt.close(fig)

cpu_count = multiprocessing.cpu_count()
pool = multiprocessing.Pool(cpu_count)
pool.starmap(processStep, [(data[i], i, shapex, shapey, type[0]) for i in [f for f in data.files if 'step' in f]])