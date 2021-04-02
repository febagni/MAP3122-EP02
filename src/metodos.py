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
import matplotlib.pyplot as plt

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
    for k in range (1,n+1):
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
    x_explicit = []
    x_explicit.append(np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t))
    x_explicit.append(np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t))
    x_explicit.append(-np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t))
    x_explicit.append(-np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t))
    return np.array(x_explicit)

def calc_error(rk4values, t0, tf, n):
    h = (tf-t0)/n
    t = t0
    errIter = []
    for k in range (1, n+1):
        #errIter = max(abs(rk4values[k] - explicit_solution(t)))
        errIter.append(max([abs(explicit_solution(t)[i] - rk4values[k][i]) for i in range (len(rk4values[k]))]))
        t = t0 + h*k
    return errIter

def get_time_array(t0,tf,n):
    time_array = []
    h = (tf-t0)/n
    t = t0
    for i in range(1,n+1):
        time_array.append(t)
        t += h
    return time_array

def plot_graph(x, y, title):
    plt.scatter(x, y, label='erros', color='r')
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title(title)
    plt.legend()
    plt.show()

    
##################################################
#exercicio1.py

x0 = [1,1,1,-1]
A = np.array([[-2,-1,-1,-2],[1,-2,2,-1],[-1,-2,-2,-1],[2,-1,1,-2]])

n = [20, 40, 80, 160, 320, 640]
xf = []
R = []

for i in range (len(n)-1):
    time_array = get_time_array(0,2,n[i])
    resultado = rk4system(x0,A,0,2,n[i])
    resultado_aux = rk4system(x0,A,0,2,n[i + 1])
    erro = calc_error(resultado,0,2,n[i])
    erro_aux = calc_error(resultado_aux,0,2,n[i + 1])
    xf.append(resultado[n[i]])
    R.append(max(erro)/max(erro_aux))
    plot_graph(time_array, erro, str("Gráfico de E(t) para n =" + str(n[i])))
    if i == len(n) - 2 : plot_graph(get_time_array(0,2,n[i+1]), erro_aux, str("Gráfico de E(t) para n =" + str(n[i+1])))

print(R)