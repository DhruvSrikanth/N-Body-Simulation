# Environment variables
PYTHON = python3
GCC = g++-11
MPI = mpiexec

# Libraries to include and link
lib_include = /opt/homebrew/Cellar/open-mpi/4.1.3/include 

lib_link_1 = /opt/homebrew/Cellar/open-mpi/4.1.3/lib
lib_link_2 = /opt/homebrew/opt/libevent/lib

# Checks
VISUALIZE = "true"

# Plot type
plot_type = "point"

# Input to cpp
n = "1000"
dt = "0.01"
NT = "1000"

# Constants
G = "6.6743e-11"

# Initialization method for the bodies
initialization = "elipsoid"

# Number of threads
num_threads = "8"

# Number of nodes
num_nodes = "4"

# Run sequence
TARGETS = output simulation visualize clean_runtime clean_text clean_png

# Run the following
target: $(TARGETS)

output:
	mkdir output

# Compile the following
nbody.out: ./nbody.cpp
	$(GCC) ./nbody.cpp -o ./nbody.out -I $(lib_include) -O3 -ffast-math -march=native -mtune=native -lm -L $(lib_link_1) -L $(lib_link_2) -lmpi -fopenmp


# Execute the following
simulation: ./nbody.out
	@echo "Simulating...\n"
	$(MPI) --bind-to none -np $(num_nodes) ./nbody.out $(n) $(dt) $(NT) $(G) $(num_threads) $(initialization)
	@echo "\nSimulation Complete!\n"

# Visualize simulation
visualize:
ifeq ($(VISUALIZE), "true")
	@echo "Visualizing...\n"
	$(PYTHON) ./visualize.py $(plot_type)
	ffmpeg -i ./output/movie.avi -strict -2 ./output/movie.mp4         
	rm ./output/movie.avi
	@echo "\nVisualization Complete!\n"	
endif

# Visualize (testing)
visualize_test:
	PYTHON visualize.py $(plot_type)
	ffmpeg -i ./output/movie.avi -strict -2 ./output/movie.mp4         
	rm ./output/movie.avi

# Clean the following
clean_runtime:
	@echo "Cleaning runtime files...\n"
	rm -rf __pycache__
	rm *.out

clean_text:
	@echo "Cleaning text files...\n"
	rm ./output/*.txt

clean_png:
	@echo "Cleaning png files...\n"
	rm ./output/*.png
	
clean_output:
	@echo "Cleaning output files...\n"
	rm -rf ./output