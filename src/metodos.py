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

# Brief : Essa funçao aplica o algoritmo Runge Kutta em uma unica iteracao
def rk4iter(x, t, h, f) :
    K1 = h*f(t,x)
    K2 = h*f(t+h/2, x + K1/2)
    K3 = h*f(t+h/2, x + K2/2)
    K4 = h*f(t+h, x + K3)
    return x + (K1 + 2*K2 + 2*K3 + K4)/6

# Brief : Essa funçao aplica o algoritmo Runge Kutta resolvendo um SPVI
def rk4system(x0, A, t0, tf, n) :
    newx = np.zeros(len(A))
    x = x0.copy()
    h = (tf-t0)/n
    t = t0
    for k in range (1,n+1) :
        for i in range (len(A)) :
            coef = A[i][i]
            sum = 0
            for j in range (len(A)) :
                if i != j : 
                    sum += A[i][j]*x[j]
            f = lambda t, x : coef*x + sum
            newx[i] = rk4iter(x[i],t,h,f)
        x = newx.copy()
    return x

##################################################
#exercicio1.py

x0 = [1,1,1,-1]
I = [0, 2]
n = 1
A = np.array([[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]])

resultado = rk4system(x0,A,0,2,1000)
print(resultado)


