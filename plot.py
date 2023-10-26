import numpy as np
import matplotlib.pyplot as plt

BIN_NUM = 30
f = open("dataBinder.txt", "r")
lines = f.readlines()
tempset = []
susset = []
f.close()


for i in range(1,len(lines), BIN_NUM+1):
    tempset.append(lines[i])
    temporary = []
    for line in lines[i+1:i+BIN_NUM+1]:
        temporary.append(eval(line))
    susset.append(np.mean(temporary))

plt.scatter(tempset,susset, marker="o", color = "IndianRed")
plt.xlabel("Temparature(T)")
plt.ylabel("Magnetization")
plt.show()