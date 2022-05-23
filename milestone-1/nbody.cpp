#include <iostream>
using namespace std;
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <omp.h>

// Data structure to hold a point in 3D Cartesian space
struct cartesian3D {
  double x;
  double y;
  double z;

};

void print_cartesian3D(cartesian3D point) {
  cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << endl;
}

// Data structure to represent each body in the simulation
struct Body {
    cartesian3D r; // position
    cartesian3D v; // velocity
    double m; // mass
};

void print_body(Body b) {
   cout << "r: (" << b.r.x << " ," << b.r.y << " ," << b.r.z << ")" << "  v:" << b.v.x << " ," << b.v.y << " ," << b.v.z << ") m: " << b.m << endl;
}

void write_to_file(Body* bodies, string filename, int n, int t) { 
    // Allocate memory for the file
    ofstream file;
    file.open(filename);

    for (int i = 0; i < n; i++) {
        if (i == 0) {
            file << "[[" << bodies[i].r.x << "," << bodies[i].r.y << "," << bodies[i].r.z << "],";
        }
        else if (i == n - 1) {
            file << "[" << bodies[i].r.x << "," << bodies[i].r.y << "," << bodies[i].r.z << "]]";
        }
        else {
            file << "[" << bodies[i].r.x << "," << bodies[i].r.y << "," << bodies[i].r.z << "],";
        }
    }

    // Release the memory for the file
    file.close();
}

void initialize_bodies(Body* bodies, int n) {
    // Set the random seed
    // srand(1); 
    int i;

    #pragma omp parallel for default(none) private(i) shared(bodies, n) schedule(guided)
    for (i = 0; i < n; i++) {
        // Initialize the bodies to a uniform distribution in the unit cube [0,1]^3
        bodies[i].r.x = (double)rand() / (double) RAND_MAX;
        bodies[i].r.y = (double)rand() / (double) RAND_MAX;
        bodies[i].r.z = (double)rand() / (double) RAND_MAX;

        // Initialize the velocity of each body to a random value in the unit cube [-1,1]^3
        bodies[i].v.x = (double)rand() / (double) RAND_MAX * 2 - 1;
        bodies[i].v.y = (double)rand() / (double) RAND_MAX * 2 - 1;
        bodies[i].v.z = (double)rand() / (double) RAND_MAX * 2 - 1;

        // Initialize the mass of each body to a random value
        bodies[i].m = 1.0;
    }    
}

void nbody(int n, double dt, int N, double G){
    // Initialize the bodies
    Body* bodies = new Body[n];
    initialize_bodies(bodies, n);
    cartesian3D F_ij;

    for (int t = 0; t < N; t++) {
        // Write the bodies to a file
        if (t % 10 == 0) {
            char filename[100];
            int dummy_var = sprintf(filename, "./output/t_%d.txt", t);
            write_to_file(bodies, filename, n, t);
        }        

        // Timing variables
        double t1;
        double t2;

        // Start the timer
        t1 = omp_get_wtime();

        // Compute the forces
        int i = 0;
        int j = 0;
        double r_mag;
        cartesian3D F_i;
        #pragma omp parallel for default(none) private(i, j, F_i, r_mag) shared(bodies, G, dt, n) schedule(guided) 
        for (i = 0; i < n; i++) {
            // Compute the force on body i
            F_i = {0, 0, 0};

            #pragma omp parallel for default(none) private(i, j, F_ij, F_i, r_mag) shared(bodies, G, n) schedule(guided)  
            for (j = 0; j < n; j++) {
                r_mag = sqrt(pow(bodies[j].r.x - bodies[i].r.x, 2) + pow(bodies[j].r.y - bodies[i].r.y, 2) + pow(bodies[j].r.z - bodies[i].r.z, 2));

                if (r_mag > 0.0) {
                    // Compute the force on body j
                    F_ij.x = G * bodies[i].m * bodies[j].m * (bodies[j].r.x - bodies[i].r.x) / pow(r_mag, 3);
                    F_ij.y = G * bodies[i].m * bodies[j].m * (bodies[j].r.y - bodies[i].r.y) / pow(r_mag, 3);
                    F_ij.z = G * bodies[i].m * bodies[j].m * (bodies[j].r.z - bodies[i].r.z) / pow(r_mag, 3);

                    // Add the force to the total force on body i
                    F_i.x += F_ij.x;
                    F_i.y += F_ij.y;
                    F_i.z += F_ij.z;
                }
            }
            // Update the velocity of body i
            bodies[i].v.x += dt * F_i.x / bodies[i].m;
            bodies[i].v.y += dt * F_i.y / bodies[i].m;
            bodies[i].v.z += dt * F_i.z / bodies[i].m;
        }
        // Update the position of each body
        #pragma omp parallel for default(none) private(i) shared(bodies, dt, n) schedule(guided) 
        for (i = 0; i < n; i++) {
            bodies[i].r.x += dt * bodies[i].v.x;
            bodies[i].r.y += dt * bodies[i].v.y;
            bodies[i].r.z += dt * bodies[i].v.z;
        }

        // Stop the timer
        t2 = omp_get_wtime();

        // Print out timestep and grind rate
        int grind_rate = floor(1.0 / (t2 - t1));
        printf("Timestep: %d - Grind Rate: %d iter/secs\n", t, grind_rate);
    }
}

int main(int argc, char** argv) {
    // Initialize variables
    int n = stoi(argv[1]); // Number of particles
    double dt = stod(argv[2]); // Timestep
    int N = stoi(argv[3]); // Number of timesteps
    double G = stod(argv[4]); // Gravitational constant
    int num_threads = stoi(argv[5]); // Number of threads

    cout << "Simulation Parameters:" << endl;
    cout << "Number of particles: " << n << endl;
    cout << "Timestep: " << dt << endl;
    cout << "Number of timesteps: " << N << endl;
    cout << "Gravitational constant: " << G << endl;

    // Set number of threads
    cout << "Number of threads = " << num_threads << "\n" << endl;
    omp_set_num_threads(num_threads);

    // Timing variables
    double ts;
    double tf;

    // Perform simulation
    ts = omp_get_wtime();
    nbody(n, dt, N, G);
    tf = omp_get_wtime();

    // Print time taken
    cout << "Time taken = " << tf-ts << " seconds" << "\n" << endl;

    return 0;

}