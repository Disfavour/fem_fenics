import numpy as np

print([-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)])
gammas = list(np.arange(-1, 4.1, 0.5))
kappas = list(np.arange(0, 2.01, 0.25))
print(gammas)
print(kappas)


alfas = [-1, 0, 1]
gammas = [-0.5, 0.2, 0, 123]

params = [(alfa, gamma)
          for alfa in alfas
          for gamma in gammas
          if not (alfa == 0 and gamma != 0)]
print(params)

alfas = [-10**i for i in range(4, -1, -1)] + [-0.5, 0.0, 0.5] + [10**i for i in range(5)]
print(alfas)
print([round(0.5 + 0.01 * 2**i, 2) for i in range(8)])