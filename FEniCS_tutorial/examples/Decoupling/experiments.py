from multiprocessing import Pool
import time

def f(x, y):
    print("hi f")
    time.sleep(5)
    return x*y

def g(x):
    print("hi g")
    time.sleep(5)
    return x*x

if __name__ == '__main__':
    with Pool(5) as p:
        #print(p.starmap(f, [(1, 1), (2, 2), (3, 3)]))
        res = p.map_async(g, [1, 2, 3])
        res1 = p.map_async(g, [4, 5, 6])
        print(res.get())
        print(res1.get())

        # res = [p.apply_async(f, t) for t in [(1, 1), (2, 2), (3, 3)]]
        # for i in res:
        #     print(i.get())

