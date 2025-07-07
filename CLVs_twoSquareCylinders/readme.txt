SC_CLV_orig.m uses Gram-Schmidt vectors and the upper-triangular R matrices previouly stored by the solver QR_code to calculate CLVs.

The resulting CLVs (CLV1,CLV2,CLV3,CLV4,CLV5 and CLV6) are stored in their respective text files.

SC_MATLAB.sh is a bash script that runs SC_CLV_orig.m on an HPC cluster.

The folder SPOD_CLV1_CLV2 contains MATLAB scripts that perform SPOD of CLV_1 and CLV_2, the two unstable CLVs of the flow. 

