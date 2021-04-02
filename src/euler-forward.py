import numpy as np

def euler_simples(y0, z0, n, I, f, g):
    x = []
    y = [y0]
    z = [z0]
    h = (I[1] - I[0])/n
    for l in range(n+1):
        x.append(I[0] + h * l)
    for l in range(1, n+1):
        K1y = f(x[l-1], y[l-1], z[l-1])
        K1z = g(x[l-1], y[l-1], z[l-1])
  
        y.append(y[l-1] + (K1y)*h)
        z.append(z[l-1] + (K1z)*h)
    Y = []
    for i in range(len(y)):
        tmp = [y[i], z[i]]
        Y.append(tmp)
        
    return np.array(Y)
    
# euler modificado
def euler_sistema(y0, z0, n, I, f, g):
    """
    y0: valor inicial de y
    z0: valor inicial de z
    n: quantos pontos teremos no dominio I
    I: lista do tipo [x0, xf]
    f: função de callback do tipo f(x, y, z)
    g: função de callback do tipo g(x, y, z)
    """
    x = []
    y = [y0]
    z = [z0]
    h = (I[1] - I[0])/n
    for l in range(n+1):
        x.append(I[0] + h * l)
    for l in range(1, n+1):
        K1y = f(x[l-1], y[l-1], z[l-1])
        K1z = g(x[l-1], y[l-1], z[l-1])
        
        K2y = f(x[l-1]+h/2, y[l-1] + h/2 * K1y, z[l-1] + h/2 * K1z)
        K2z = g(x[l-1]+h/2, y[l-1] + h/2 * K1y, z[l-1] + h/2 * K1z)
        
        y.append(y[l-1] + (K2y)/2)
        z.append(z[l-1] + (K2z)/2)
    Y = []
    for i in range(len(y)):
        tmp = [y[i], z[i]]
        Y.append(tmp)
        
    return np.array(Y)

y0 = 4
z0 = 2
I = [0, 1]
n = 2
f = lambda x, y, z: z
g = lambda x, y, z: -x*z-2*y
resultado = euler_simples(y0, z0, n, I, f, g)
print(resultado)