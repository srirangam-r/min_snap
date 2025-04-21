import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# load JSON
with open("min_snap_2d_coeffs.json") as f:
    j = json.load(f)

times    = np.array(j["times"])
coeffs_x = np.array(j["coeffs_x"]).reshape(-1,8)
coeffs_y = np.array(j["coeffs_y"]).reshape(-1,8)

# sampling
t_final = times[-1]
dt      = 0.01
t_vals  = np.arange(0, t_final+dt, dt)
x_vals  = np.zeros_like(t_vals)
y_vals  = np.zeros_like(t_vals)

# evaluate piecewise 
for i in range(len(coeffs_x)):
    if i < len(coeffs_x)-1:
        mask = (t_vals >= times[i]) & (t_vals <  times[i+1])
    else:
        mask = (t_vals >= times[i]) & (t_vals <= times[i+1])
    tau = t_vals[mask] - times[i]
    for j in range(8):
        x_vals[mask] += coeffs_x[i,j] * tau**j
        y_vals[mask] += coeffs_y[i,j] * tau**j

fig, ax = plt.subplots()
ax.set_aspect('equal','box')
ax.set_xlim(x_vals.min(), x_vals.max())
ax.set_ylim(y_vals.min(), y_vals.max())
ax.plot(x_vals, y_vals, color='lightgray', lw=1)

point, = ax.plot([], [], 'ro')
trail, = ax.plot([], [], 'b-', lw=2)

def init():
    # initialize point and trail
    point.set_data([], [])
    trail.set_data([], [])
    return point, trail

def update(frame):
    # update point and trail
    point.set_data((x_vals[frame],), (y_vals[frame],))
    trail.set_data(x_vals[:frame+1], y_vals[:frame+1])
    return point, trail

anim = FuncAnimation(
    fig, update,
    frames=len(t_vals),
    init_func=init,
    blit=False,
    interval=20
)

# save 
writer = PillowWriter(fps=30)
anim.save("min_snap_2d.gif", writer=writer)
print("Saved min_snap_2d.gif")
