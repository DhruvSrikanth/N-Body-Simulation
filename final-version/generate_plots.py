import numpy as np
import matplotlib.pyplot as plt

def generate_strong_scaling_plot(n, t1, t2, filename):
    """
    Generates a plot of the strong scaling results.
    """

    n = np.array(n)
    t1 = np.array(t1)
    t2 = np.array(t2)

    # Plot the data
    plt.plot(n, t1, '-o', label='Pure MPI')
    plt.plot(n, t2, '-o', label='Hybrid MPI+OpenMP')
    # Add legened
    plt.legend(loc="upper left")
    # Set the labels
    plt.xlabel("n (number of nodes)")
    plt.ylabel("T(n) - Runtime (secs)")
    # Set the title
    plt.title("Strong Scaling Analysis")
    # Save the plot
    plt.savefig(filename)
    plt.clf()

n1 = [1,2,4,8]
time1 = [657.2141, 328.7498, 164.802, 81.9527]
time2 = [27.4820, 13.8322, 6.9990, 3.5570]
generate_strong_scaling_plot(n1, time1, time2, "strong_scaling_comparsion.png")