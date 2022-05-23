# Project 4 - N Body Problem

## Local System Information:
### Processor:

```
Chip:               Apple M1 Pro
Number of Cores:	8 (6 performance and 2 efficiency)
```

### Compiler:

```
Compiler:           GNU C++ compiler - g++
Version:            11
```

# Milestone 1 - Serial CPU Implementation (with OpenMP)

The purpose of this portion of the project is to develop a serial, cpu-based algorithm for simulating the n-body problem. 

## Run Instructions:
To run the n-body algorithm written in `C++`, follow these steps.

1. Enter the correct directory - 
```
cd milestone-1
```

2. Specify the input `n` **(number of bodies)**, `dt` **(timestep)**, `N` **(number of timesteps)** and `num_threads` **(number of OpenMP threads)** in the Makefile.

3. Build and run - 
```
make
```

This will generate the executable, run the code for the inputs specified in the Makefile, and generate the plot for the render.

4. Clean the temporary files - 
```
make clean_text
make clean_png
```

The outputs can be found in the `milestone-1/output` directory as `movie.mp4`.

## Results:

### Inputs - 
```
1. n (number of bodies) = 102,400
2. dt (timestep) = 0.1
3. N (number of timesteps) = 500
4. num_threads (number of OpenMP threads) = 8
```

### Time Taken - 

```
Time Taken = 10.3677 seconds
```

### Inputs - 
```
1. n (number of bodies) = 1000
2. dt (timestep) = 0.001
3. N (number of timesteps) = 1000
4. num_threads (number of OpenMP threads) = 8
```

### Results - 

The animation of the simulation with the above input parameters can be found [here](https://drive.google.com/file/d/1a8nYa6tsxLdRzUjwXkKhkyGafxZbqthp/view?usp=sharing) under `movie.mp4`.


