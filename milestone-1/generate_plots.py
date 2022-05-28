import numpy as np
import matplotlib.pyplot as plt

def generate_strong_scaling_plot(n, time, analysis_num):
    """
    Generates a plot of the strong scaling results.
    """
    t1 = time[0]
    Sn = [t1/t for t in time]

    n = np.array(n)
    Sn = np.array(Sn)

    # Plot the data
    plt.plot(n, Sn, '-o')
    # Set the labels
    plt.xlabel("n (number of cores)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Strong Scaling Analysis")
    # Save the plot
    plt.savefig("./strong_scaling_{}.png".format(analysis_num))
    plt.clf()

n1 = [1,2,4,6,8]
time1 = [17.3392, 8.9305, 4.5864, 3.1169, 2.8383]
generate_strong_scaling_plot(n1, time1, "1")