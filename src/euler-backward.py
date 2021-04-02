import matplotlib.pyplot as plt
import numpy as np

def partial_derivative(t, u, partial_in, f, step):
    step_u = u.copy()
    step_u[partial_in] += step
    return (f(t,step_u) - f(t,u))/step

def newton_iter(t, u, f, h, last_u):
    n = len(u)
    jacobian = np.identity(n)
    G = np.zeros(len(u))
    new_t = t + h
    for i in range(n):
        for j in range(n):
            jacobian[i][j] -= h*partial_derivative(new_t, u, j, f[i], 0.001)
    inv_jacobian = np.linalg.inv(jacobian)
    for i in range(len(u)):
        G[i] = u[i] - h * f[i](t,u) - last_u[i]
    return u - np.matmul(inv_jacobian, G) # == new_u

def implicit_euler_iter(u, t, h, f):
    newton_u = u.copy()
    euler_u = u.copy()
    for _ in range(7) :
        newton_u = newton_iter(t, newton_u, f, h, u)
    #print("u e newton u")
    #print(u,newton_u)
    for i in range(len(u)):
        euler_u[i] = u + h * f[i](t, newton_u)
    return euler_u

def implicit_euler_system(u, f, t0, tf, n):
    u_values = []
    new_u = u.copy()
    u_values.append(new_u)
    h = (tf-t0)/n
    t = t0
    for _ in range(1, n+1):
        new_u = implicit_euler_iter(new_u,t,h,f)
        u_values.append(new_u)
        t += h
    return u_values
    
u = np.array([-8.79])
f = []
f.append(lambda t,x : 2*t + (x-t*t)*(x-t*t))

resposta_euler_implicito = implicit_euler_system(u, f, 1.1, 3.0, 5000)
print(resposta_euler_implicito[5000])