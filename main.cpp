#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip> // For formatting the table

// Constants
const size_t NUM_PARTICLES = 1000000;
const double HBAR = 1.0545718e-34; // Reduced Planck's constant

// Array of Structs (AoS) implementation
struct ParticleAoS {
    double position[3];  // x, y, z
    double momentum[3];  // px, py, pz
    double spin[2];      // spin components (simplified)
};

// Struct of Arrays (SoA) implementation
struct ParticleSoA {
    std::vector<double> position_x;
    std::vector<double> position_y;
    std::vector<double> position_z;
    std::vector<double> momentum_x;
    std::vector<double> momentum_y;
    std::vector<double> momentum_z;
    std::vector<double> spin_x;
    std::vector<double> spin_y;

    ParticleSoA(size_t n) {
        position_x.resize(n); position_y.resize(n); position_z.resize(n);
        momentum_x.resize(n); momentum_y.resize(n); momentum_z.resize(n);
        spin_x.resize(n); spin_y.resize(n);
    }
};

// Function to initialize particles with random values
void initializeAoS(std::vector<ParticleAoS>& particles) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (auto& p : particles) {
        for (int i = 0; i < 3; i++) {
            p.position[i] = dis(gen);
            p.momentum[i] = dis(gen);
        }
        p.spin[0] = dis(gen);
        p.spin[1] = dis(gen);
    }
}

void initializeSoA(ParticleSoA& particles) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (size_t i = 0; i < NUM_PARTICLES; i++) {
        particles.position_x[i] = dis(gen);
        particles.position_y[i] = dis(gen);
        particles.position_z[i] = dis(gen);
        particles.momentum_x[i] = dis(gen);
        particles.momentum_y[i] = dis(gen);
        particles.momentum_z[i] = dis(gen);
        particles.spin_x[i] = dis(gen);
        particles.spin_y[i] = dis(gen);
    }
}

// Computation 1: Calculate total phase (position * momentum / hbar)
double computePhaseAoS(const std::vector<ParticleAoS>& particles) {
    double total_phase = 0.0;
    for (const auto& p : particles) {
        for (int i = 0; i < 3; i++) {
            total_phase += p.position[i] * p.momentum[i] / HBAR;
        }
    }
    return total_phase;
}

double computePhaseSoA(const ParticleSoA& particles) {
    double total_phase = 0.0;
    for (size_t i = 0; i < NUM_PARTICLES; i++) {
        total_phase += (particles.position_x[i] * particles.momentum_x[i] +
                       particles.position_y[i] * particles.momentum_y[i] +
                       particles.position_z[i] * particles.momentum_z[i]) / HBAR;
    }
    return total_phase;
}

// Computation 2: Update spin states (simple rotation)
void updateSpinAoS(std::vector<ParticleAoS>& particles) {
    for (auto& p : particles) {
        double temp = p.spin[0];
        p.spin[0] = p.spin[1];
        p.spin[1] = -temp;
    }
}

void updateSpinSoA(ParticleSoA& particles) {
    for (size_t i = 0; i < NUM_PARTICLES; i++) {
        double temp = particles.spin_x[i];
        particles.spin_x[i] = particles.spin_y[i];
        particles.spin_y[i] = -temp;
    }
}

int main() {
    // Memory usage estimation
    size_t aos_size = NUM_PARTICLES * sizeof(ParticleAoS);
    size_t soa_size = NUM_PARTICLES * 8 * sizeof(double); // 8 arrays
    double aos_mem_mb = aos_size / 1024.0 / 1024.0;
    double soa_mem_mb = soa_size / 1024.0 / 1024.0;

    // Initialize AoS
    std::vector<ParticleAoS> particles_aos(NUM_PARTICLES);
    initializeAoS(particles_aos);

    // Initialize SoA
    ParticleSoA particles_soa(NUM_PARTICLES);
    initializeSoA(particles_soa);

    // Benchmark AoS
    auto start = std::chrono::high_resolution_clock::now();
    double phase_aos = computePhaseAoS(particles_aos);
    updateSpinAoS(particles_aos);
    auto end = std::chrono::high_resolution_clock::now();
    auto aos_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // Benchmark SoA
    start = std::chrono::high_resolution_clock::now();
    double phase_soa = computePhaseSoA(particles_soa);
    updateSpinSoA(particles_soa);
    end = std::chrono::high_resolution_clock::now();
    auto soa_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // Output results in the requested table format
    std::cout << "\nNums of particles = " << NUM_PARTICLES << "\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "                " << std::setw(15) << "AoS" << std::setw(15) << "SoA" << "\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "Memory usage    " << std::fixed << std::setprecision(2) << std::setw(15) << aos_mem_mb 
              << std::setw(15) << soa_mem_mb << "\n";
    std::cout << "Time            " << std::setw(15) << aos_time 
              << std::setw(15) << soa_time << "\n";
    std::cout << "Phase           " << std::scientific << std::setprecision(2) << std::setw(15) << phase_aos 
              << std::setw(15) << phase_soa << "\n";
    std::cout << "-----------------------------------------------------------------\n";

    return 0;
}