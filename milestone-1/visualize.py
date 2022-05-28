import os
import sys

import ast

import numpy as np
import matplotlib.pyplot as plt
import cv2

# plt.style.use('seaborn-poster')

def read_file(filename):
    '''
    Read the file and return the data.
    '''
    with open(filename, "r") as f:
        file_str = f.read()
        data = ast.literal_eval(file_str)
        grid = np.array(data)
    return grid

def plot_point(point, ax):
    '''
    Plot the point.
    '''
    x, y, z = point
    ax.scatter3D(x, y, z, cmap='winter')
    return ax

def plot_surface(grid, ax):
    grid = np.array(grid)
    X = grid[:,0]
    Y = grid[:,1]
    Z = grid[:,2]

    ax.plot_trisurf(X, Y, Z, cmap='winter', edgecolor='none')
    
    return ax
    
def generate_movie(plot_flag="surface"):
    '''
    Plot and save the data.
    '''
    step = 10
    fig = plt.figure(figsize = (10,10))

    NT = len(os.listdir('./output/')) * step
    
    for i in range(0, NT, step):
        try:
            filename = './output/t_' + str(i) + '.txt'
            data = read_file(filename)
            print("Plotting Simulation Sample : " + str(i))

            ax = fig.add_subplot(111, projection='3d')
            ax.grid(True)
            # ax.axis('off')
            # ax.set_facecolor('xkcd:black')

            if plot_flag == "surface":
                ax = plot_surface(data, ax)
                ax.view_init(30, 60)
            elif plot_flag == "point":
                for point in data:
                    ax = plot_point(point, ax)
                ax.view_init(27, 32)
            else:
                raise ValueError("Invalid plot flag")
            
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            ax.set_title('Simulation Sample : ' + str(i))

            plt.savefig('./output/t_' + str(i) + '.png')
            plt.clf()

        except FileNotFoundError:
            continue

    movie_name = "./output/movie.avi"
    files = ['./output/t_' + str(i) + '.png' for i in range(0, NT, step)]
    frames = [cv2.imread(file) for file in files]
    
    # save the frames as a movie
    frame = frames[0]
    height, width, _ = frame.shape  
    video = cv2.VideoWriter(movie_name, 0, 1, (width, height)) 
    for frame in frames: 
        video.write(frame) 
      
    cv2.destroyAllWindows() 
    video.release()
        
        

if __name__ == '__main__':
    gen_mov = True

    # Read in from the command line
    plot_flag = sys.argv[1]

    if gen_mov:
        generate_movie(plot_flag=plot_flag)

