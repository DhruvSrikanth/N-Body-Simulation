import numpy as np
import matplotlib.pyplot as plt
import ast
import os
import cv2


plt.style.use('seaborn-poster')

def read_file(filename):
    '''
    Read the file and return the data.
    '''
    with open(filename, "r") as f:
        file_str = f.read()
        data = ast.literal_eval(file_str)
        grid = np.array(data)
    return grid

def generate_movie():
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

            for j in range(0, len(data)):
                x, y, z = data[j]
                ax.scatter3D(x, y, z)
            
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_zlim(0, 1)

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

    if gen_mov:
        generate_movie()

