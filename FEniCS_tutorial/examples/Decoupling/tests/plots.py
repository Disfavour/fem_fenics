import numpy as np
import matplotlib.pyplot as plt




nonlinear_p1 = np.load("P1/data/nonlinear.npy")
jac1_p1 = np.load("P1/data/jac1.npy")
jac2_p1 = np.load("P1/data/jac2.npy")
jac5_p1 = np.load("P1/data/jac5.npy")
seidel1_p1 = np.load("P1/data/seidel1.npy")
seidel2_p1 = np.load("P1/data/seidel2.npy")
seidel5_p1 = np.load("P1/data/seidel5.npy")

nonlinear_p2 = np.load("P2/data/nonlinear.npy")
jac1_p2 = np.load("P2/data/jac1.npy")
jac2_p2 = np.load("P2/data/jac2.npy")
jac5_p2 = np.load("P2/data/jac5.npy")
seidel1_p2 = np.load("P2/data/seidel1.npy")
seidel2_p2 = np.load("P2/data/seidel2.npy")
seidel5_p2 = np.load("P2/data/seidel5.npy")

time = np.linspace(0, 5, nonlinear_p1.size)

plt.figure(dpi=100, figsize=(6.4, 4), tight_layout=True)


def jac_p1():
    plt.plot(time, nonlinear_p1, label="nonlinear_p1")
    plt.plot(time, jac1_p1, label="jac1_p1")
    plt.plot(time, jac2_p1, label="jac2_p1")
    plt.plot(time, jac5_p1, label="jac5_p1")


def seidel_p1():
    plt.plot(time, nonlinear_p1, label="nonlinear_p1")
    plt.plot(time, seidel1_p1, label="seidel1_p1")
    plt.plot(time, seidel2_p1, label="seidel2_p1")
    plt.plot(time, seidel5_p1, label="seidel5_p1")


def jac_p2():
    plt.plot(time, nonlinear_p2, label="nonlinear_p2")
    plt.plot(time, jac1_p2, label="jac1_p2")
    plt.plot(time, jac2_p2, label="jac2_p2")
    plt.plot(time, jac5_p2, label="jac5_p2")


def seidel_p2():
    plt.plot(time, nonlinear_p2, label="nonlinear_p2")
    plt.plot(time, seidel1_p2, label="seidel1_p2")
    plt.plot(time, seidel2_p2, label="seidel2_p2")
    plt.plot(time, seidel5_p2, label="seidel5_p2")


# plt.plot(time, jac5_p1 - nonlinear_p1, label="jac5_p1 - nonlinear_p1")
# plt.plot(time, seidel5_p1 - nonlinear_p1, label="seidel5_p1 - nonlinear_p1")

#plt.plot(time, jac5_p2 - nonlinear_p2, label="jac5_p2 - nonlinear_p2")
#plt.plot(time, seidel5_p2 - nonlinear_p2, label="seidel5_p2 - nonlinear_p2")

jac_p1()

plt.legend()
plt.grid()
plt.xlabel("Time")
plt.ylabel("Energy")

plt.savefig("dif_p2.pdf")

plt.show()
