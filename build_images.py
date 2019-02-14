import numpy as np
import sys
import matplotlib.pyplot as plt
def cell_to_pos(msx, msy, pos):
    return (pos % msx, int(pos / msx))

filename = sys.argv[1]

nb_step  = sys.argv[2]
displ = len(nb_step)
nb_step = int(nb_step)

data = np.load(filename)

shapex, shapey = data['shape']

print(shapex, shapey)

grid = np.zeros((shapex, shapey))
for i in range(nb_step):
    image = data['step-%d' % i]
    assert image.shape[0] == shapey*shapex
    fig, ax = plt.subplots(1, 1, figsize=(20,20))
    for gid in image:
        x, y = cell_to_pos(shapex, shapey, gid)
        grid[x][y] = 1
    ax.imshow(grid)

    fig.savefig(("{0:0>%d}_water_dummy.jpg" % displ).format(i))
    plt.close(fig)





