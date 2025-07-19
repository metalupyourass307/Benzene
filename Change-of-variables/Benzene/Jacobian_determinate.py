import numpy as np
from x_bar_related import *
from y_bar_related import *
from Phi import *

lambda_z = 0.01
def density_grid(x, y, x_alpha, y_alpha, lambda_z):
    return np.sum(1 / np.sqrt(lambda_z**2 + (x - x_alpha)**2 + (y - y_alpha)**2 )) 

def d_Phi(t):
    logterm=np.log((1+t)/(1-t))
    return 2/((1-t**2)(np.sqrt(1+logterm**2)))

def J_det(v, w, x, y, lambda_x, x_h, x_l, y_h, y_l, x_alpha, y_alpha, lambda_y,):
    a_x = alpha_bar_x(x_alpha, lambda_x, x_h, x_l)
    alpha_y_u = alpha_y(x, y_h, y_l, x_alpha, y_alpha, lambda_y)
    density_g = density_grid(x, y, x_alpha, y_alpha, lambda_z)
    
    return (a_x * alpha_y_u * Phi(v) * Phi(w))/density_g