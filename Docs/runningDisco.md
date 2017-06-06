** Running Disco **

========

**** First Time Setup ****

Disco relies on the HDF5 (for efficient i/o) and MPI (for parallism) libraries. The location of these libraries, as well as other system-specifc options, are specified by `Makefile_dir.in`. We provide a template in `Makefile_dir.in.template`. First copy the template into a new file:

    $ cp Makefile_dir.in.template Makefile_dir.in

Then modify the entries as necessary. You should not need to change this file much (if at all) once you have Disco running.
