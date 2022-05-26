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

void initialize_bodies(Body* bodies, int n, string initialization_type) {
    // Set the random seed
    // srand(1); 
    int i;
    double phi, theta;

    #pragma omp parallel for default(none) private(i, phi, theta) shared(bodies, n, initialization_type) schedule(guided)
    for (i = 0; i < n; i++) {
        // sample random phi and theta
        phi = 2 * M_PI * (double)rand() / (double)RAND_MAX;
        theta = acos((double)rand() / (double)RAND_MAX);

        if (initialization_type == "random") {
            // set position, velocity, and mass to random values
            bodies[i].r.x = (double)rand() / (double)RAND_MAX;
            bodies[i].r.y = (double)rand() / (double)RAND_MAX;
            bodies[i].r.z = (double)rand() / (double)RAND_MAX;

            bodies[i].v.x = (double)rand() / (double)RAND_MAX;
            bodies[i].v.y = (double)rand() / (double)RAND_MAX;
            bodies[i].v.z = (double)rand() / (double)RAND_MAX;

            bodies[i].m = (double)rand() / (double)RAND_MAX;
        }
        else if (initialization_type == "uniform") {
            // set position, velocity, and mass to uniform values
            bodies[i].r.x = (double)(i + 1) / (double)n;
            bodies[i].r.y = (double)(i + 1) / (double)n;
            bodies[i].r.z = (double)(i + 1) / (double)n;

            bodies[i].v.x = (double)(i + 1) / (double)n;
            bodies[i].v.y = (double)(i + 1) / (double)n;
            bodies[i].v.z = (double)(i + 1) / (double)n;

            bodies[i].m = (double)(i + 1) / (double)n;
        }
        else if (initialization_type == "elipsoid") {
            // Set the position to follow parametrized elipse
            bodies[i].r.x = cos(theta) * sin(phi);
            bodies[i].r.y = sin(theta) * sin(phi);
            bodies[i].r.z = cos(phi);

            // Set the velocity to follow parametrized elipse
            bodies[i].v.x = -sin(theta) * sin(phi);
            bodies[i].v.y = cos(theta) * sin(phi);
            bodies[i].v.z = -sin(phi);

            // Set the mass to be a function of the radius to the center of the elipse
            bodies[i].m = 1.0 / (1.0 + bodies[i].r.x * bodies[i].r.x + bodies[i].r.y * bodies[i].r.y + bodies[i].r.z * bodies[i].r.z);
        }
        else if (initialization_type == "galaxy") {
            // Set the position to follow parametrized spiral
            bodies[i].r.x = cos(phi);
            bodies[i].r.y = (i + 1)*sin(phi);
            bodies[i].r.z = 1;

            // Set the velocity to follow parametrized spiral
            bodies[i].v.x = -sin(phi);
            bodies[i].v.y = cos(phi);
            bodies[i].v.z = 1.0;

            // Set the mass to be a function of the radius to the center of the spiral
            bodies[i].m = 1.0 / (1.0 + bodies[i].r.x * bodies[i].r.x + bodies[i].r.y * bodies[i].r.y + bodies[i].r.z * bodies[i].r.z);
        }
        else {
            assert("Invalid initialization type");
            exit(1);
        }        
    }    
}

void nbody(int n, double dt, int N, double G, string initialization_type){
    // Initialize the bodies
    Body* bodies = new Body[n];
    initialize_bodies(bodies, n, initialization_type);
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
        #pragma omp parallel for default(none) private(i, j, F_ij, F_i, r_mag) shared(bodies, G, dt, n) schedule(guided) 
        for (i = 0; i < n; i++) {
            // Compute the force on body i
            F_i = {0, 0, 0};
            
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
    string initialization_type = argv[6]; // Initialization type

    cout << "Simulation Parameters:" << endl;
    cout << "Number of particles: " << n << endl;
    cout << "Timestep: " << dt << endl;
    cout << "Number of timesteps: " << N << endl;
    cout << "Gravitational constant: " << G << endl;
    cout << "Initialization type: " << initialization_type << endl;

    // Set number of threads
    omp_set_num_threads(num_threads);
    cout << "Number of threads = " << num_threads << "\n" << endl;

    // Timing variables
    double ts;
    double tf;

    // Perform simulation
    ts = omp_get_wtime();
    nbody(n, dt, N, G, initialization_type);
    tf = omp_get_wtime();

    // Print time taken
    cout << "Time taken = " << tf-ts << " seconds" << "\n" << endl;

    return 0;

}