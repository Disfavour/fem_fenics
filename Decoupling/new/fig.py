import matplotlib.pyplot as plt
import numpy as np


plt.figure(figsize=(6.4, 3.6), tight_layout=True)

old = np.load(f"data/old.npy")
new = np.load(f"data/new.npy")

plt.plot(old[:, 0], old[:, 1], color='blue', linestyle="solid", label=r"old")
plt.plot(old[:, 0], new[:, 1], color='red', linestyle="solid", label=r"new")

plt.legend()
plt.grid()
plt.xlim(old[:, 0].min(), old[:, 0].max())
plt.xlabel(r"$t$")
plt.ylabel(r"$m$")

for i in [".pdf", ".png"]:
    plt.savefig("plots/m" + i, dpi=300, bbox_inches='tight')

plt.figure(figsize=(6.4, 3.6), tight_layout=True)

plt.plot(old[:, 0], old[:, 2], color='blue', linestyle="solid", label=r"old")
plt.plot(old[:, 0], new[:, 2], color='red', linestyle="solid", label=r"new")

plt.legend()
plt.grid()
plt.xlim(old[:, 0].min(), old[:, 0].max())
plt.xlabel(r"$t$")
plt.ylabel(r"$E$")

for i in [".pdf", ".png"]:
    plt.savefig("plots/E" + i, dpi=300, bbox_inches='tight')
