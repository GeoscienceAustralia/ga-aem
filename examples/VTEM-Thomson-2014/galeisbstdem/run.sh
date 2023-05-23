#!/bin/tcsh

# Note: The control files galeisbstdem_ascii_dfn.con and galeisbstdem_ascii_column.con have the same inversion settings. However the _dfn.con version uses the ASEGGDF2 .dfn file variable names and the _column.con version uses column numbers to define fields.
# e.g. In galeisbstdem_ascii_column.con TX_Height = Column 30
#      In galeisbstdem_ascii_dfn.con    TX_Height = emloop_height

# Use a single standalone CPU
galeisbstdem.exe galeisbstdem_ascii_dfn.con
# galeisbstdem.exe galeisbstdem_ascii_column.con

# Use 4 MPI processes
# mpirun -np 4 galeisbstdem.exe galeisbstdem_ascii_dfn.con
# mpirun -np 4 galeisbstdem.exe galeisbstdem_ascii_column.con

# Use 4 OpenMP threads
# galeisbstdem.exe galeisbstdem_ascii_dfn.con 4
# galeisbstdem.exe galeisbstdem_ascii_column.con 4

