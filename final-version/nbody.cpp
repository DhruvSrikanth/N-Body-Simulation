#include <iostream>
using namespace std;
#include <stdio.h>
#include <fstream>
#include <assert.h>
#include <string>
#include <math.h>
#include <omp.h>
#include <mpi.h>

// Dimension of node arrangement in the MPI grid
#define DIMENSION 1

// Data structure to hold a point in 3D Cartesian space
struct cartesian3D {
  double x;
  double y;
  double z;

};

void print_cartesian3D(cartesian3D point) {
  cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << endl;
}

cartesian3D add_vectors(cartesian3D a, cartesian3D b) {
    cartesian3D c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;

    return c;
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

int local_to_global_1d_offset(int mype, int nprocs) {
    int offset = mype * nprocs;
    return offset;
}

void copy_array(Body *source, Body *destination, int n) {
    int i;
    #pragma omp parallel for default(none) private(i) shared(source, destination, n) schedule(guided)  
    for (i = 0; i < n; i++) {
        destination[i] = source[i];
    }
}

void initialize_bodies(Body* bodies, int n_local, string initialization_type, int mype, int nprocs) {
    // Set the random seed
    // srand(1); 
    int i;
    int global_i;
    int offset = local_to_global_1d_offset(mype, nprocs);

    double phi, theta;

    #pragma omp parallel for default(none) private(i, global_i, phi, theta) shared(bodies, n_local, initialization_type, offset) schedule(guided)
    for (i = 0; i < n_local; i++) {
        // Compute the global index of the body for some types of initialization
        global_i = i + offset;

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
            bodies[i].r.x = (double)(global_i + 1) / (double)n_local;
            bodies[i].r.y = (double)(global_i + 1) / (double)n_local;
            bodies[i].r.z = (double)(global_i + 1) / (double)n_local;

            bodies[i].v.x = (double)(global_i + 1) / (double)n_local;
            bodies[i].v.y = (double)(global_i + 1) / (double)n_local;
            bodies[i].v.z = (double)(global_i + 1) / (double)n_local;

            bodies[i].m = (double)(global_i + 1) / (double)n_local;
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
            bodies[i].r.y = (global_i + 1)*sin(phi);
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

double distance(Body b1, Body b2) {
    double res;
    double dx = b1.r.x - b2.r.x;
    double dy = b1.r.y - b2.r.y;
    double dz = b1.r.z - b2.r.z;
    res = sqrt(dx * dx + dy * dy + dz * dz);
    return res;
}

cartesian3D compute_force(Body b1, Body b2, double distance_ij, double G) {
    // Compute the force on body 1 due to body 2
    cartesian3D force;

    force.x = G * b1.m * b2.m * (b2.r.x - b1.r.x) / pow(distance_ij, 3);
    force.y = G * b1.m * b2.m * (b2.r.y - b1.r.y) / pow(distance_ij, 3);
    force.z = G * b1.m * b2.m * (b2.r.z - b1.r.z) / pow(distance_ij, 3);

    return force;
}

Body update_velocity(Body b, cartesian3D F, double dt) {
    // Compute the velocity of body b
    b.v.x += dt * F.x / b.m;
    b.v.y += dt * F.y / b.m;
    b.v.z += dt * F.z / b.m;
    return b;
}

void update_positions(Body *bodies, int n_local, double dt) {
    int i;
    #pragma omp parallel for default(none) private(i) shared(bodies, dt, n_local) schedule(guided) 
    for (i = 0; i < n_local; i++) {
        bodies[i].r.x += dt * bodies[i].v.x;
        bodies[i].r.y += dt * bodies[i].v.y;
        bodies[i].r.z += dt * bodies[i].v.z;
    }
}

void nbody(int n_local, double dt, int N, double G, string initialization_type, int mype, int nprocs, MPI_Comm comm1d, MPI_Datatype MPI_Body){
    // Initialize the bodies
    Body* bodies_local = new Body[n_local];
    Body* bodies_remote = new Body[n_local];
    Body* bodies_buffer = new Body[n_local];

    initialize_bodies(bodies_local, n_local, initialization_type, mype, nprocs);

    // Send local to remote and receive remote to local
    for (int r = 0; r < nprocs; r++) {
        if (mype == 0) {
            initialize_bodies(bodies_remote, n_local, initialization_type, r, nprocs);
            // Send to remote
            MPI_Send(bodies_remote, n_local, MPI_Body, r, 99, comm1d);
        }
        else {
            // Receive from remote
            MPI_Recv(bodies_local, n_local, MPI_Body, 0, 99, comm1d, MPI_STATUS_IGNORE);
        }
    }

    cartesian3D F_ij;

    for (int t = 0; t < N; t++) {
        // Write the bodies to a file
        if (t % 10 == 0) {
            char filename[100];
            int dummy_var = sprintf(filename, "./output/t_%d.txt", t);
            // write_to_file(bodies, filename, n_local, t);
        }

        // Timing variables
        double t1;
        double t2;

        // Start the timer
        t1 = omp_get_wtime();

        // Copy local bodies to remote
        copy_array(bodies_local, bodies_remote, n_local);

        // Pipeline for computing forces on local bodies from other local and remote bodies
        for (int r = 0; r < nprocs; r++) {
            // Compute the forces
            int i = 0;
            int j = 0;
            double distance_ij;
            cartesian3D F_i;

            // For each local body
            #pragma omp parallel for default(none) private(i, j, F_ij, F_i, distance_ij) shared(bodies_local, bodies_remote, G, dt, n_local) schedule(guided) 
            for (i = 0; i < n_local; i++) {
                // Compute the force on body i
                F_i = {0, 0, 0};

                // For each remote body
                for (j = 0; j < n_local; j++) {
                    distance_ij = distance(bodies_local[j], bodies_remote[i]);

                    if (distance_ij > 0.0) {
                        F_ij = compute_force(bodies_local[i], bodies_remote[j], distance_ij, G);
                        // Add the force to the total force on body i
                        F_i = add_vectors(F_i, F_ij);
                    }
                }
                
                // Update the velocity of body i
                bodies_local[i] = update_velocity(bodies_local[i], F_i, dt);
            }
            // Send Bremote buffer to left neighbor rank
            MPI_Status status;
            MPI_Sendrecv(bodies_remote, n_local, MPI_Body, r - 1, 99, bodies_buffer, n_local, MPI_Body, r + 1, MPI_ANY_TAG, comm1d, &status);
            // Receive new Bremote buffer from right neighbor rank
            copy_array(bodies_buffer, bodies_remote, n_local);
        }

        // Update the position of each body
        update_positions(bodies_local, n_local, dt);

        // Stop the timer
        t2 = omp_get_wtime();

        // Print out timestep and grind rate
        int grind_rate = floor(1.0 / (t2 - t1));
        printf("Timestep: %d - Grind Rate: %d iter/secs\n", t, grind_rate);
    }
}

int main(int argc, char** argv) {
    // MPI initialization
    int mype, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    // Create a type for struct cartesian 3D
    int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype MPI_Cartesian3D;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(cartesian3D, x);
    offsets[1] = offsetof(cartesian3D, y);
    offsets[2] = offsetof(cartesian3D, z);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_Cartesian3D);
    MPI_Type_commit(&MPI_Cartesian3D);

    // Create a type for struct body
    nitems = 3;
    blocklengths[0] = 1;
    blocklengths[1] = 1;
    blocklengths[2] = 1;
    types[0] = MPI_Cartesian3D;
    types[1] = MPI_Cartesian3D;
    types[2] = MPI_DOUBLE;
    MPI_Datatype MPI_Body;
    offsets[0] = offsetof(Body, r);
    offsets[1] = offsetof(Body, v);
    offsets[2] = offsetof(Body, m);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_Body);
    MPI_Type_commit(&MPI_Body);


    // MPI Cartesian Grid Creation (1D)
    int dims[DIMENSION], periodic[DIMENSION], coords[DIMENSION];
    MPI_Comm comm1d;
    dims[0] = nprocs;
    periodic[0] = 1;

    // Create Cartesian Communicator
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periodic, 1, &comm1d);

    // Extract this MPI rank's N-dimensional coordinates from its place in the MPI Cartesian grid
    MPI_Cart_coords(comm1d, mype, DIMENSION, coords);

    // Determine 1D neighbor ranks for this MPI rank
    int left, right;
    MPI_Cart_shift(comm1d, 0, 1, &left, &right);
    // cout << "Rank " << mype << " has left neighbor " << left << " and right neighbor " << right << endl;

    // Initialize variables
    int n_global = stoi(argv[1]); // Number of particles

    // Check whether the physical domain can be divided evenly among the MPI ranks
    assert(n_global % nprocs == 0);

    // Number of particles per MPI rank
    int n_local = n_global / nprocs;

    double dt = stod(argv[2]); // Timestep
    int N = stoi(argv[3]); // Number of timesteps
    double G = stod(argv[4]); // Gravitational constant
    int num_threads = stoi(argv[5]); // Number of threads
    string initialization_type = argv[6]; // Initialization type

    if (mype == 0) {
        cout << "Simulation Parameters:" << endl;
        cout << "Number of particles: " << n_global << endl;
        cout << "Timestep: " << dt << endl;
        cout << "Number of timesteps: " << N << endl;
        cout << "Gravitational constant: " << G << endl;
        cout << "Initialization type: " << initialization_type << endl;

        // Set number of threads
        omp_set_num_threads(num_threads);
        cout << "Number of threads = " << num_threads << "\n" << endl;
    }

    // Timing variables
    double ts;
    double tf;

    // Perform simulation
    ts = omp_get_wtime();
    nbody(n_local, dt, N, G, initialization_type, mype, nprocs, comm1d, MPI_Body);
    tf = omp_get_wtime();

    // Print time taken
    if (mype == 0) {
        cout << "Time taken = " << tf-ts << " seconds" << "\n" << endl;
    }

    MPI_Type_free(&MPI_Body);
    MPI_Type_free(&MPI_Cartesian3D);

    // MPI finalization
    MPI_Finalize();

    return 0;

}