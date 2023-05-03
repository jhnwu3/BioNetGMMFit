import numpy as np

xDir = "X/"
yDir = "Y/"
xFile = "Xslim.csv"
yFiles = ["Yt0slim.csv","Yt1slim.csv", "Yt2slim.csv", "Yt3slim.csv", "Yt4slim.csv"]

x = np.genfromtxt(xDir + xFile, delimiter=',')
xInt = x.astype(int)
np.savetxt(xDir + 'X.csv', xInt, delimiter=',', fmt='%d') 

for file in yFiles:
    y = np.genfromtxt(yDir + file,delimiter=',')
    yInt = y.astype(int)
    np.savetxt(yDir + file[:-8] + ".csv", yInt, delimiter=',', fmt='%d')