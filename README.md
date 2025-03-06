# Comparing Structure of Arrays (SoA) vs Array of Structures (AoS)

## Problem Description

In large-scale quantum physics simulations, tracking the properties of thousands or even millions of particles efficiently is crucial. Each particle has several physical properties, such as position, momentum, and spin. These properties can be stored in memory using different approaches.

The **Array of Structures (AoS)** method stores all properties of each particle together. This makes accessing all properties of a single particle easy, but accessing the same property across many particles (for example, all x-coordinates) requires skipping over other properties, making memory access less efficient.

The **Structure of Arrays (SoA)** method takes a different approach: it stores each property in its own separate array. This layout is more efficient when performing calculations that involve the same property across many particles, as the data is stored contiguously in memory.

In large-scale simulations, choosing the right layout has a direct impact on performance, affecting both memory bandwidth usage and cache efficiency.

---

## Example Output

```
Nums of particles = 1000000
-----------------------------------------------------------------
                            AoS            SoA     Difference
-----------------------------------------------------------------
Memory usage              61.04          61.04           0.00
Time (ms)                    19             15              4
Phase                 -5.20e+36      -1.06e+36
-----------------------------------------------------------------
```

---

## Why SoA is Faster

The **Structure of Arrays (SoA)** approach is faster primarily because of how modern processors access memory. Processors read memory in blocks (cache lines), and operations that require sequential access to large amounts of data are much faster when that data is stored contiguously.

In the SoA layout, all values for a single property (like position or spin) are stored in one long, continuous block of memory. This allows the processor to efficiently load large portions into cache, reducing the number of memory accesses and improving overall performance. It also enables better use of **SIMD (Single Instruction, Multiple Data)** hardware, which can process multiple values at once.

In contrast, the **Array of Structures (AoS)** layout scatters different properties together, so accessing all particles' x-coordinates requires skipping over unrelated properties (y, z, momentum, spin, etc.). This leads to inefficient cache use, more cache misses, and slower performance.

In simulations where millions of particles are updated every timestep, this difference adds up, making SoA a significantly better choice for memory-intensive numerical simulations.

---

## How to Build and Run

### Clone the Repository

```bash
git clone https://github.com/LyudmilaKostanyan/AoS-vs-SoA.git
cd AoS-vs-SoA
```

### Build the Project

```bash
cmake -S . -B build
cmake --build build
```

### Run the Program

Navigate into the `build` directory:

```bash
cd build
```

Run the executable with the number of particles as a parameter:

```bash
./main.exe --num_particles 10000
```

Replace `10000` with the desired number of particles.
