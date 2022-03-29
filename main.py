import numpy as np
import matplotlib.pyplot as plt

with open("E:\CompMathLab2Python\P.txt", 'r') as inp:
    lst = inp.read().split()
    lst = [[float(n) for n in x.split()] for x in lst]
print(lst)
x = len(lst)
j = 0.15/x

t = np.arange(-10, 10, 0.2)

plt.title('Co = 0.005, t=0.015, E')
plt.plot(t, lst)
plt.savefig("E.png")
plt.show()
plt.clf()
