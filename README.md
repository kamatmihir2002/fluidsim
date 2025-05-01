# multiphase_fluid_sim

A multiphase fluid simulation with parallel Schwarz domain decomposition with coarse correction for my parallel systems (CSC 548) project at North Carolina State University, Spring 2025.

# Building

To build this project, run:
```
make -B all
```

To build this in debug mode, run:
```
make -B alld
```

### NOTE: Please don't forget the -B flag, otherwise make gives an error.

### NOTE 2: As of now, only the debug configuration generates the visualization csv files. Visualization generation is disabled in the normal configuration for timing measurement.

# Running

An example command to run a fluid sim:
```
./fluidsim -prefix <result file prefix> -grid-extents 2.0 -grid-subdiv 128 -parallel-count 4 -steps 100 -num_particles 0
```

This performs the fluid simulation on a grid spanning from 0.0 to 2.0 both in x and y directions, subdivided into 128 elements and split into 4 * 4 parallel domains, run for 100 timesteps.

# Visualization
You can visualize the results using the ``` fluid_visualizer.py ``` python script. You need numpy, pandas and matplotlib to run the visualizer. Once these are installed, you can run this as 
```
./python3 fluid_visualizer.py <fluid grid file>.csv <None>/<particle file>.csv
```