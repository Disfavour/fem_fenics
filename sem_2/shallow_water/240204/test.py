import numpy as np

print([-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)])
