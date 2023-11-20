import numpy as np
import matplotlib.pyplot as plt

BIN_NUM = 30

def write(filename,dotcolor):
    f = open(".\\"+filename, "r")
    lines = f.readlines()
    tempset = []
    susset = []
    sigmaset = []
    f.close()


    for i in range(1,len(lines), BIN_NUM+1):
        tempset.append(lines[i])
        temporary = []
        for line in lines[i+1:i+BIN_NUM+1]:
            temporary.append(eval(line))
        susset.append(np.mean(temporary))
        sigmaset.append(np.var(temporary) / (BIN_NUM-1))
        print(np.var(temporary,ddof=1)) 

    plt.errorbar(tempset,susset, yerr=sigmaset, marker="o", color = dotcolor)

if __name__ == "__main__":
    # write("sW12.txt","IndianRed")
    write("sL12.txt","Green")
    write("sL16.txt","RoyalBlue")
    write("sL20.txt","Purple")
    write("sL24.txt","Orange")
    # write("sinMetro16.txt","Pink")
    plt.xlabel("Temparature(T)")
    plt.ylabel("Binder(L=10)")
    plt.legend(["L12","L16","L20","L24","W24","W28","M16"])
    # plt.legend(["LB Thermolize and Metro measure","pure LB measurement","Limited size LB","Wolff","PureMetropolis","PureMetro2"])
    plt.show()




