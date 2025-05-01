import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as anim

df = pd.read_csv(sys.argv[1], header=None)
Lp = np.max(df[0])
print(Lp)

L = []
for tstep in range(Lp):
    timestep_data = df[df[0] == tstep].drop(df.columns[0], axis=1)
    sortedp = timestep_data.sort_values(by=timestep_data.columns[0])
    print(sortedp)
    num_particles = len(timestep_data[timestep_data.columns[0]])
    positions = []
    positions.append(sortedp[sortedp.columns[1]].astype(np.float32))
    positions.append(sortedp[sortedp.columns[2]].astype(np.float32))
    

    K = np.array(positions)
    L.append(K)

L = np.array(L)
print(L.shape)

fig, ax = plt.subplots(1, 1)

ax.set_xlim((0.0, 1.0))
ax.set_ylim((0.0, 1.0))
A,  = ax.plot(np.zeros(10), np.zeros(10), marker='o', linestyle='none')

oldc = None
def animate(i): 
    
    A.set_xdata(i[0])
    A.set_ydata(i[1])
    return A,
    

animp = anim.FuncAnimation(fig, animate, frames=L, interval=50)
fig.tight_layout()
plt.show()



