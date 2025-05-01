import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import pandas as pd

import sys


df = pd.read_csv(sys.argv[1], header=None)
Lp = np.max(df[0])

if (len(sys.argv) > 2):
    df2 = pd.read_csv(sys.argv[2], header=None)
    Lp2 = np.max(df2[0])

pmax = -100000000
print(Lp)
L = []
L2 = []
for tstep in range(Lp):
    timestep_data = df[df[0] == tstep].drop(df.columns[0], axis=1)
    xval = timestep_data[timestep_data[1] == 'x']
    yval = timestep_data[timestep_data[1] == 'y']
    p = timestep_data[timestep_data[1] == 'p']
    positions = []
    if (len(sys.argv) > 2):
        timestep_data2 = df2[df2[0] == tstep].drop(df2.columns[0], axis=1)
        sortedp = timestep_data2.sort_values(by=timestep_data2.columns[0])
        num_particles = len(timestep_data[timestep_data2.columns[0]])
        positions.append(sortedp[sortedp.columns[1]].astype(np.float32))
        positions.append(sortedp[sortedp.columns[2]].astype(np.float32))

    K = [np.array(xval.drop(xval.columns[0], axis=1).values, dtype=np.float64),
            np.array(yval.drop(yval.columns[0], axis=1).values, dtype=np.float64),
            np.array(p.drop(yval.columns[0], axis=1).values, dtype=np.float64)]
    
    if (len(sys.argv) > 2):
        K.append(sortedp[sortedp.columns[1]].astype(np.float32))
        K.append(sortedp[sortedp.columns[2]].astype(np.float32))

    fpmax = np.max(K[2])
    if fpmax > pmax:
        pmax = fpmax
    
    
    L.append(K)

a = L[0][0].shape[0]
b = L[0][0].shape[1]
fig, ax = plt.subplots(1, 1)
A = None
if (len(sys.argv) > 2):
    A,  = ax.plot(np.zeros(10), np.zeros(10), marker='o', linestyle='none')

Q = ax.quiver(np.ones((a, b)), np.ones((a, b)), scale=50.0, cmap="viridis")


# V = ax.contour([])
oldc = None
def animate(i): 
    global A
    global oldc
    if oldc is not None:
        oldc.remove()
    Q.set_UVC(i[0], i[1])
    oldc = ax.contourf(i[2] / pmax, alpha=0.5, levels=100, antialiased=True, cmap="inferno")
    
    if (len(sys.argv) > 2):
        A.set_xdata(i[3] * (a - 1))
        A.set_ydata(i[4] * (b - 1))
    else:
        A = None
    return Q, oldc, A

animp = anim.FuncAnimation(fig, animate, frames=L, interval=10)
fig.tight_layout()
plt.show()
