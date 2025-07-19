import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from x_bar_related import *
from y_bar_related import *
from z_bar_related import *
from mpl_toolkits.mplot3d import Axes3D

# load settings from txt file
df = pd.read_csv('/Users/laiyuxuan/Desktop/Change-of-variables/Single Atom/Single_Atom.txt', delim_whitespace=True, skiprows=2, header=None)
df.columns = ['Element', 'X', 'Y', 'Z']

# atom coordinates in Cartesian space
x_alpha = df['X'].values
lambda_x = 0.01
y_alpha = df['Y'].values
lambda_y = 0.01
z_alpha = df['Z'].values
lambda_z = 0.01

# Lines in 2D Cartesian space, each has a distance 10 a.u. to the nearest atom
x_h=5
x_l=-5
y_h=5
y_l=-5
z_h=5
z_l=-5

# generate 64*64 uniform grid in 2D space spanned by the v and w coordinates
u = np.array([(i - 32) / 32 for i in range(64)])
v = np.array([(j - 32) / 32 for j in range(64)])
w = np.array([(k - 32) / 32 for k in range(64)])

# using Newton-Raphson method to numerically solve x
def newton_raphson_x(u, x0=0.0, tol=1e-10, max_iter=50):
    x = x0
    target = x_bar(u, x_alpha, lambda_x, x_h, x_l)
    for _ in range(max_iter):
        fx = x_bar_x(x, x_alpha, lambda_x) - target
        dfx = dx_bar(x, x_alpha, lambda_x)
        dx = -fx / dfx
        x += dx
        if abs(dx) < tol:
            break
    return x

# using Newton-Raphson method to numerically solve y
def newton_raphson_y(v, x, y0=0.0, tol=1e-10, max_iter=50):
    target = y_bar(v, x, y_h, y_l, x_alpha, y_alpha, lambda_y)
    y = y0
    for _ in range(max_iter):
        fy = y_bar_xy(x, y, x_alpha, y_alpha, lambda_y) - target
        dfy = dy_bar(x, y, y_alpha, x_alpha, lambda_y)
        dy = -fy / dfy
        y += dy
        if abs(dy) < tol:
            break
    return y

def newton_raphson_z(w, x, y, z0=0.0, tol=1e-8, max_iter=50):
    target = z_bar(w, x, y, z_h, z_l, x_alpha, y_alpha, z_alpha, lambda_x, lambda_y, lambda_z)
    z = z0
    for _ in range(max_iter):
        fz = z_bar_xyz(x, y, z, x_alpha, y_alpha, z_alpha, lambda_x, lambda_y, lambda_z) - target
        dfz = dz_bar(x, y, z, x_alpha, y_alpha, z_alpha, lambda_x, lambda_y, lambda_z)
        dz = -fz / dfz
        z += dz
        if abs(dz) < tol:
            break
    return z

x_vals = []
y_vals = []
z_vals = []
u_vals = np.array([(i - 32) / 32 for i in range(64)])
v_vals = np.array([(j - 32) / 32 for j in range(64)])
w_vals = np.array([(j - 32) / 32 for j in range(64)])

'''
basic idea for Newton-Raphson method here: 
1. calculate all the x_bar_targets and y_bar_targets from genrated 64*64 (v,w), via x_bar and y_bar defined in x_bar_related and y_bar_related
2. x_bar_x and y_bar_xy shows the connection between (x,y)and (x_bar,y_bar). Now that x_bar_targets and y_bar_targets has gained, we can get
    x_targets and y_targets by Newton-Raphson method, with the help of x_bar_x, y_bar_xy and their derivatives dx_bar and dy_bar.
'''

for i, u in enumerate(u_vals):
    x_sol = newton_raphson_x(u)
    for j, v in enumerate(v_vals):
        y_sol = newton_raphson_y(v, x_sol)
        for k, w in enumerate(w_vals):
            z_sol = newton_raphson_z(w, x_sol, y_sol)
            print(f"u = {v:.3f}, v = {w:.3f}, w = {w:.3f} -> x = {x_sol:.6f}, y = {y_sol:.6f}, z = {z_sol:.6f}")
            if -10 <= x_sol <= 10 and -7.5 <= y_sol <= 10:
                x_vals.append(x_sol)
                y_vals.append(y_sol)
                z_vals.append(z_sol)


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')  

ax.scatter(x_vals, y_vals, z_vals, s=2, c=z_vals, cmap='viridis')  
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Mapped (x, y, z) from (u, v, w)')

plt.tight_layout()
plt.show()