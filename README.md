# vasp-mdgrabber
This program grabs a position, force configuration from the VASP OUTPUT file, either, an abinitio molecular dynamics OUTCAR or any sort of VASP calculation. It is written in FORTRAN and can be further extended to grab such data from other DFT input for further processing. 

- First of all compile, the file `processoutcar.f90` using `gfortran processoutcar.f90 -o your_output_name`, where `your_output_name` is the name of objective file you want to name

- To run you have to have, the OUTCAR in the same location of the binary or `your_output_name`. If you have this simply run `./your_output_name` from your working directory
