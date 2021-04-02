##############################################################
# File      :   metodos.py
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
def rk4system(x0, A, t0, tf, n):
    newx = np.zeros(len(A))
    x = x0.copy()
    rk4values = []
    rk4values.append(x)
    h = (tf-t0)/n
    t = t0
    for _ in range (1,n+1):
        for i in range (len(A)):
            coef = A[i][i]
            sum = 0
            for j in range (len(A)):
                if i != j : 
                    sum += A[i][j]*x[j]
            f = lambda t, x : coef*x + sum
            newx[i] = rk4iter(x[i],t,h,f)
        t = t + h
        x = newx.copy()
        rk4values.append(x)
    return np.array(rk4values)
    
# Inserir comparacao entre solucao explicita e atual iteracao = obter erro
        # Comparar com o maior erro (preciso do max Err)
        # Armazenar qual for maior
        #E = max([abs(x_explicit[i] - x[i]) for i in range (4)])

def explicit_solution(t):
    x_explicit = np.array([
                 np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
                 np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t),
                 -np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
                 -np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t)])
    return x_explicit

def calc_error(rk4values, t0, tf, n):
    h = (tf-t0)/n
    t = t0
    errIter = []
    for k in range (1, n+1):
        #errIter = max(abs(rk4values[k] - explicit_solution(t)))
        errIter.append(max([abs(explicit_solution(t)[i] - rk4values[k][i]) for i in range (len(rk4values[k]))]))
        t = t0 + h*k
    return errIter