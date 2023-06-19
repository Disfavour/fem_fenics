import numpy as np
import matplotlib.pyplot as plt


nonlinear = np.load("nonlinear.npy")
decoupled = np.load("decoupled.npy")
decoupled_nabla = np.load("decoupled_nabla.npy")
decoupled_jac = np.load("decoupled_jac.npy")
decoupled_seidel = np.load("decoupled_seidel.npy")

decoupled2 = np.load("decoupled2.npy")
decoupled_jac2 = np.load("decoupled_jac2.npy")
decoupled_seidel2 = np.load("decoupled_seidel2.npy")

plt.plot(nonlinear, label="nonlinear")

# plt.plot(decoupled, label="decoupled")
# #plt.plot(decoupled_nabla, label="decoupled_nabla")
# plt.plot(decoupled_jac, label="decoupled_jac")
# plt.plot(decoupled_seidel, label="decoupled_seidel")

#plt.plot(decoupled2, label="decoupled2")
plt.plot(decoupled_jac2, label="decoupled_jac2")
plt.plot(decoupled_seidel2, label="decoupled_seidel2")

plt.legend()
plt.show()
