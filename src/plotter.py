import matplotlib.pyplot as plt

def get_time_array(t0,tf,n):
    time_array = []
    h = (tf-t0)/n
    t = t0
    for _ in range(1,n+1):
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