##############################################################
# File      :   metodos.vhd
# Project   :   MAP3122 - EP02 - Metodos numericos para 
#                                resolucao de EDOs
# Date      :   April/2021
##############################################################
# Author    :   Felipe Bagni
# Email     :   febagni@usp.br
# NUSP      :   11257571
##############################################################
# Author    :   Gabriel Yugo Kishida
# Email     :   gabriel.kishida@usp.br
# NUSP      :   11257647
##############################################################

import numpy as np 

def runge_kutta_4(x0, n, I, f):
    """
    x0: valor inicial de x
    n: quantos pontos teremos no dominio I
    I: lista do tipo [x0, xf]
    f: função do tipo f(t, x) callback que com base em valor de t e x, retorna resultado 
    """
    t = []
    x = []
    x.append(x0)
    newx = []
    h = (I[1] - I[0])/n
    
    for l in range(n+1):
        t.append(I[0] + h * l)

    for l in range(n):
        for j in range(4):
            print("l,j ",l,j)
            #K1 = f(t[l], x[l], j)
            #K2 = f(t[l] + h/2, x[l] + K1*h/2, j)
            #K3 = f(t[l] + h/2, x[l] + K2*h/2, j)
            #K4 = f(t[l] + h, x[l] + K3*h, j)
            #newx.append(x[l][j] + (K1+2*K2+2*K3+K4)*h/6)

            K1 = f(t[l], x[l], j)
            x[l][j] += K1*h/2
            K2 = f(t[l] + h/2, x[l], j)
            x[l][j] += K2*h/2
            K3 = f(t[l] + h/2, x[l], j)
            x[l][j] += K3*h
            K4 = f(t[l] + h, x[l], j)
            x[l][j] -= K1*h/2 + K2*h/2 + K3*h
            newx.append(x[l][j] + (K1+2*K2+2*K3+K4)*h/6)
            print(newx)
        x.append(newx) 
        newx.clear()
    return np.array(x)

##################################################
#exercicio1.py

x0 = [1,1,1,-1]
I = [0, 2]
n = 1
A = np.array([[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]])

print(np.matmul(A[0],x0))

def f(t, x, i):
    A = np.array([[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]]) 
    return np.matmul(A[i],x)

#f = []
#for i in range(len(A)) :
#    f.append(lambda t, x: np.matmul(A[i],x))
    
resultado = runge_kutta_4(x0, n, I, f)
print(resultado)
