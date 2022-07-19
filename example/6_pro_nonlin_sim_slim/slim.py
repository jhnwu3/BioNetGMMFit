import numpy as np
def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)

my_data = np.genfromtxt('Y.csv', delimiter=',')
# my_data1 = np.genfromtxt('Yt1.csv', delimiter=',')
# my_data2 = np.genfromtxt('Yt2.csv', delimiter=',')
# my_data3 = np.genfromtxt('Yt3.csv', delimiter=',')
# my_data4 = np.genfromtxt('Yt4.csv', delimiter=',')
my_data = my_data[:5000,:]
# my_data1 = my_data1[:5000,:]
# my_data2 = my_data2[:5000,:]
# my_data3 = my_data3[:5000,:]
# my_data4 = my_data4[:5000,:]

# my_data = trunc(my_data,decs=4)
# print(my_data)
np.savetxt('Y.csv',my_data, delimiter=',',fmt='%1.4f')
# np.savetxt('Yt0slim.csv',my_data, delimiter=',',fmt='%1.4f')
# np.savetxt('Yt1slim.csv',my_data1, delimiter=',',fmt='%1.4f')
# np.savetxt('Yt2slim.csv',my_data2, delimiter=',',fmt='%1.4f')
# np.savetxt('Yt3slim.csv',my_data3, delimiter=',',fmt='%1.4f')
# np.savetxt('Yt4slim.csv',my_data4, delimiter=',',fmt='%1.4f')