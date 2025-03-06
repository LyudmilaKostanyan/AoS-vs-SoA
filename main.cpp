#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <string>
#include "kaizen.h"

size_t num_particles = 1000000;
const double HBAR = 1.0545718e-34;

struct ParticleAoS {
    double position[3];
    double momentum[3];
    double spin[2];
};

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

    for (size_t i = 0; i < num_particles; i++) {
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
    for (size_t i = 0; i < num_particles; i++) {
        total_phase += (particles.position_x[i] * particles.momentum_x[i] +
                       particles.position_y[i] * particles.momentum_y[i] +
                       particles.position_z[i] * particles.momentum_z[i]) / HBAR;
    }
    return total_phase;
}

void updateSpinAoS(std::vector<ParticleAoS>& particles) {
    for (auto& p : particles) {
        double temp = p.spin[0];
        p.spin[0] = p.spin[1];
        p.spin[1] = -temp;
    }
}

void updateSpinSoA(ParticleSoA& particles) {
    for (size_t i = 0; i < num_particles; i++) {
        double temp = particles.spin_x[i];
        particles.spin_x[i] = particles.spin_y[i];
        particles.spin_y[i] = -temp;
    }
}

int main(int argc, char** argv) {
    zen::cmd_args args(argv, argc);
    auto particles_options = args.get_options("--num_particles");
    if (!particles_options.size() || std::stoi(particles_options[0]) <= 0)
        std::cerr << zen::color::red("num_particles has been set to the default value of 1000000.") << std::endl
                  << zen::color::red("To set your own value, pass it as an argument to your executable:") << std::endl
                  << "--num_particles *your_value*" << std::endl
                  << zen::color::red("The value must be greater than 0.") << std::endl;
    else
        num_particles = atoi(particles_options[0].c_str());
    
    size_t aos_size = num_particles * sizeof(ParticleAoS);
    size_t soa_size = num_particles * 8 * sizeof(double);
    double aos_mem_mb = aos_size / 1024.0 / 1024.0;
    double soa_mem_mb = soa_size / 1024.0 / 1024.0;
    double mem_diff = aos_mem_mb - soa_mem_mb;

    std::vector<ParticleAoS> particles_aos(num_particles);
    initializeAoS(particles_aos);
    
    ParticleSoA particles_soa(num_particles);
    initializeSoA(particles_soa);

    zen::timer timer;
    timer.start();
    double phase_aos = computePhaseAoS(particles_aos);
    updateSpinAoS(particles_aos);
    timer.stop();
    auto aos_time = timer.duration<zen::timer::msec>().count();

    timer.start();
    double phase_soa = computePhaseSoA(particles_soa);
    updateSpinSoA(particles_soa);
    timer.stop();
    auto soa_time = timer.duration<zen::timer::msec>().count();
    auto time_diff = aos_time - soa_time;

    std::cout << "\nNums of particles = " << num_particles << "\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "                " << std::setw(15) << "AoS" << std::setw(15) << "SoA" << std::setw(15) << "Difference" << "\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "Memory usage    " << std::fixed << std::setprecision(2) << std::setw(15) << aos_mem_mb 
              << std::setw(15) << soa_mem_mb << std::setw(15) << mem_diff << "\n";
    std::cout << "Time (ms)       " << std::setw(15) << aos_time 
              << std::setw(15) << soa_time << std::setw(15) << time_diff << "\n";
    std::cout << "Phase           " << std::scientific << std::setprecision(2) << std::setw(15) << phase_aos 
              << std::setw(15) << phase_soa << std::setw(15) << "\n";
    std::cout << "-----------------------------------------------------------------\n";

    return 0;
}