import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from x_bar_related import *
from y_bar_related import *
from Phi import *

# load settings from txt file
df = pd.read_csv('Benzene.txt', delim_whitespace=True, skiprows=2, header=None)
df.columns = ['Element', 'X', 'Y', 'Z']

# atom coordinates in Cartesian space
x_alpha = df['X'].values
lambda_x = 0.01
y_alpha = df['Y'].values
lambda_y = 0.01

# Lines in 2D Cartesian space, each has a distance 10 a.u. to the nearest atom
x_h=16.0471
x_l=-13.4015
y_h=16.3772
y_l=-11.8003

# generate 64*64 uniform grid in 2D space spanned by the v and w coordinates
v = np.array([(i - 32) / 32 for i in range(64)])
w = np.array([(j - 32) / 32 for j in range(64)])

# using Newton-Raphson method to numerically solve x
def newton_raphson_x(v, x0=0.0, tol=1e-10, max_iter=50):
    x = x0
    target = x_bar(v, x_alpha, lambda_x, x_h, x_l)
    for _ in range(max_iter):
        fx = x_bar_x(x, x_alpha, lambda_x) - target
        dfx = dx_bar(x, x_alpha, lambda_x)
        dx = -fx / dfx
        x += dx
        if abs(dx) < tol:
            break
    return x

# using Newton-Raphson method to numerically solve y
def newton_raphson_y(w, x, y0=0.0, tol=1e-10, max_iter=50):
    target = y_bar(w, x, y_h, y_l, x_alpha, y_alpha, lambda_y)
    y = y0
    for _ in range(max_iter):
        fy = y_bar_xy(x, y, x_alpha, y_alpha, lambda_y) - target
        dfy = dy_bar(x, y, y_alpha, x_alpha, lambda_y)
        dy = -fy / dfy
        y += dy
        if abs(dy) < tol:
            break
    return y

x_vals = []
y_vals = []
v_vals = np.array([(i - 32) / 32 for i in range(64)])
w_vals = np.array([(j - 32) / 32 for j in range(64)])

'''
basic idea for Newton-Raphson method here: 
1. calculate all the x_bar_targets and y_bar_targets from genrated 64*64 (v,w), via x_bar and y_bar defined in x_bar_related and y_bar_related
2. x_bar_x and y_bar_xy shows the connection between (x,y)and (x_bar,y_bar). Now that x_bar_targets and y_bar_targets has gained, we can get
    x_targets and y_targets by Newton-Raphson method, with the help of x_bar_x, y_bar_xy and their derivatives dx_bar and dy_bar.
'''

for i, v in enumerate(v_vals):
    x_sol = newton_raphson_x(v)
    for j, w in enumerate(w_vals):
        y_sol = newton_raphson_y(w, x_sol)
        print(f"v = {v:.3f}, w = {w:.3f} -> x = {x_sol:.6f}, y = {y_sol:.6f}")
        if -10 <= x_sol <= 10 and -7.5 <= y_sol <= 10:
            x_vals.append(x_sol)
            y_vals.append(y_sol)

plt.figure(figsize=(6,6))
plt.scatter(x_vals, y_vals, s=2)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Mapped (x, y) Grid from (v, w)")
plt.axis("equal")
plt.grid(True)
plt.tight_layout()
plt.show()