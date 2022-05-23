# Environment variables
PYTHON = python3
GCC = g++-11

# Checks
VISUALIZE = "true"

# Input to cpp
n = "500"
dt = "0.1"
NT = "2000"
G = "6.6743e-11"

# Number of threads
num_threads = "8"

# Run sequence
TARGETS = output simulation visualize clean_runtime clean_text clean_png

# Run the following
target: $(TARGETS)

output:
	mkdir output

# Compile the following
nbody.out: ./nbody.cpp
	$(GCC) ./nbody.cpp -o ./nbody.out -O3 -ffast-math -march=native -mtune=native -lm -fopenmp

# Execute the following
simulation: ./nbody.out
	@echo "Simulating...\n"
	./nbody.out $(n) $(dt) $(NT) $(G) $(num_threads)
	@echo "\nSimulation Complete!\n"

# Visualize simulation
visualize:
ifeq ($(VISUALIZE), "true")
	@echo "Visualizing...\n"
	$(PYTHON) ./visualize.py
	ffmpeg -i ./output/movie.avi -strict -2 ./output/movie.mp4         
	rm ./output/movie.avi
	@echo "\nVisualization Complete!\n"	
endif

# Visualize (testing)
visualize_test:
	PYTHON visualize.py
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