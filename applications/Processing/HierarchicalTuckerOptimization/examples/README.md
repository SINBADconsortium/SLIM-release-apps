# Missing receiver interpolation of 3D frequency slices using Hierarchical Tucker optimization.

## Scripts in this directory
1. interp4D.m - 3D frequency slice interpolation using the Hierarchical Tucker optimization code using dense linear algebra routines. This script supports both serial, for data fitting in to one core, and parallel modes, for distributed data. Can handle real and complex data.

2. interp4Dsparse.m - same as the above, but uses sparse linear algebra routines (more efficient when the amount of data available is much smaller than say 1% of the total data volume). Serial and parallel modes available. Uses a mex file implementation (supports interpolating real data only for the moment).

## Running
From this directory's parent
1. Run the 'startup.m' script.

2. Fetch data from the FTP server by calling 'scons' in the /data directory.

3. Run the 'interp4D.m' script in the /examples directory to run the interpolation script.

4. Run the 'interp4Dview.m' script in the /docs directory to view the results.

