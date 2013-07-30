import numpy as np
import matplotlib.pyplot as plt

def main():
    x = np.arange(-1.0, 1.0, 0.01)
    for p in xrange(10):
        y = x ** p * (x - 1) * (x + 1)
        plt.plot(x, y, label = "p = %d" % (p))

    plt.legend()
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()
