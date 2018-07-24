import numpy as np
import matplotlib.pyplot as plt 
data = np.loadtxt("block_2d_g0_8-64_2")
plt.figure()
plt.plot(data)
plt.title("d = 2, g = 0, $kappa$ = 0.124843945\n $L_X = 8$, $L_T = 64$, nblock = 5000, blocks = 10000")
plt.savefig("2d_8-64_5000_10000")
#size = len(data) - 20
#corr = sum(data[20:])
#corr /= size
#data_sq = np.array(data) ** 2
#corr_sq = sum(data_sq[20:])
#corr_sq /= size
#err = np.sqrt((corr_sq - corr**2)/size)
#print(corr, err)
data2 = [sum(data[i*10:i*10+10])/10 for i in range(500)]
plt.figure()
plt.plot(data2)
plt.title("d = 2, g = 0, $kappa$ = 0.124843945\n $L_X = 8$, $L_T = 64$, nblock = 500, blocks = 100000")
plt.savefig("2d_8-64_500_100000")

data3 = [sum(data2[i*10:i*10+10])/10 for i in range(50)]
plt.figure()
plt.plot(data3)
plt.title("d = 2, g = 0, $kappa$ = 0.124843945\n $L_X = 8$, $L_T = 64$, nblock = 50, blocks = 1000000")
plt.savefig("2d_8-64_50_1000000")

data4 = [sum(data2[i*5:i*5+5])/5 for i in range(10)]
plt.figure()
plt.plot(data4)
plt.title("d = 2, g = 0, $kappa$ = 0.124843945\n $L_X = 8$, $L_T = 64$, nblock = 10, blocks = 5000000")
plt.savefig("2d_8-64_10_5000000")

