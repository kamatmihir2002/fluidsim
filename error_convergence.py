import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
df = pd.read_csv(sys.argv[1], header=None)
errors = np.array(df[df.columns[0]], dtype = np.float64)
iters = np.array(df[df.columns[1]], dtype = np.int64)
fig, ax1 = plt.subplots()
ax1.set_ylim(0.0, 0.5)

ax2 = ax1.twinx()
ax2.set_ylim(0.0, 200)
plt.title(sys.argv[2])
P1, = ax1.plot(errors, color="tab:red", label="Error convergence")
P2, = ax2.plot(iters, color="tab:green", label="Iterations")
ax1.set_ylabel("Error magnitude")
ax2.set_ylabel("Iteration count")
plt.legend(handles=[P1, P2])
plt.show()
