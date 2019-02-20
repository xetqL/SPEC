import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
from math import log
f = open('perf.log', 'r')
payload = f.read()
lines = payload.split('\n')[:-1]
data = [json.loads(l.split(']: ')[1]) for l in lines]
loads = pd.Series([d['loads'] for d in data])
slopes= pd.Series([d['slopes'] for d in data])
loads2 = loads.apply(lambda x: [ float((v-min(x))/(max(x)-min(x)))  for v in x])
slopes2 = slopes.apply(lambda x: [float((v-min(x))/(max(x)-min(x))) if (max(x) > 0) else 0 for v in x])

loadarr = np.array(loads2.values.tolist())
slopearr = np.array(slopes2.values.tolist())
loadarr1 = np.array(loads.values.tolist())

hists = loads.apply(lambda x: [float(v) for v in np.histogram(x, bins=15)[0]])
slopehists = slopes.apply(lambda x: [float(v) for v in np.histogram(x, bins=15)[0]])
[print(max(v)) for v in loadarr1]
histarr = np.array(hists.values.tolist())
slopehistarr = np.array(slopehists.values.tolist())
li = [d['load_imbalance'] for d in data]
fig, ax = plt.subplots(1, 4)
ax[0].plot(li)
ax[1].imshow(slopehistarr)
ax[2].imshow(loadarr)
ax[3].imshow(slopearr)
plt.show(fig)

