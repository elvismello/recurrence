"""
Just having sure I'm making things as they are supposed to be made

"""


import numpy as np
import matplotlib.pyplot as plt


def logistic_trajectory(x, r=4.0, iter=1000):
    """
    Creates a logistic trajectory based in the initial conditions given.
    """
    
    trajectory = []
    for i in range(iter):
        trajectory.append(x)
        x = r * x * (1 - x)
    
    return trajectory

if __name__ == '__main__':
    
    x0 = 0.4
    epsilon = 0.01
    iterations = 1000

    # generating data
    trajectory = logistic_trajectory(x0, iter=iterations)
    time = [i for i in range(iterations)]


    # creating recurrence matrix
    matrix = np.zeros((iterations, iterations))
    for i in range(iterations):
        for j in range(iterations):
            if abs(trajectory[i] -trajectory[j]) < epsilon:
                matrix[i, j] = 1
            else:
                matrix[i, j] = 0
    
    plt.imshow(matrix, cmap='Greys', origin='lower')
    plt.savefig('pythonTest.png', dpi=600)

    plt.close()
    plt.plot(time, trajectory, lw=0.5)
    plt.show()



