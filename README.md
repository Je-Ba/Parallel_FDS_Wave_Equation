## Solving Wave Equation using MPI
The repository includes a serial and parallel version of a finite difference solver which solver the two dimensional wave equation. The parallel version decomposes the domain into rectangular strips and utilizes MPI Non-Blocking communications to share data among the processes. This leads to a significant speed-up compared to the serial implementation when run on multiple cores. Other changes to the initial code include the use of arrays instead of 2-D vectors which leads to better performance as well as more boundary conditions and a easy to use user input file. A post-processing code  can generate gif of the simulation which makes is easy to visualize the results. Example outputs for different boundary conditions
are shown below:


<img align="left" src="https://user-images.githubusercontent.com/72440497/115837135-290bc000-a410-11eb-90fa-1d48835975e8.gif" width="350" height="350"/>
<img align="right" src="https://user-images.githubusercontent.com/72440497/115837630-b3ecba80-a410-11eb-8a39-d23b7755588d.gif" width="350" height="350"/>
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

<img align="left" src="https://user-images.githubusercontent.com/72440497/115837435-7ee06800-a410-11eb-8aba-e8810f439c86.gif" width="350" height="350"/>
<img align="right" src="https://user-images.githubusercontent.com/72440497/115837446-8142c200-a410-11eb-93e2-af1df637b353.gif" width="350" height="350"/>
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />





# Install
Simply clone this repository on to your local and compile the c++ code in Parallel_Wave_solver.cpp using Microsoft Visual Studios with MPI enabled or on a HPC system compile with: <br />
`mpic++ Parallel_Wave_solver.cpp -lm` <br />
To use the post-processing, please use Python 3 and pip install all packages in the requirements.txt file. <br />

# Use
After the code is complied, simply run the executable but make sure that the User_Input.txt file is in the same folder. The user can change all simulation inputs such as size of the domain, end time, boundary condition or the location of a Neumann Barrier in this text file which will be read in during the set-up phase of the code. This makes it very easy for the user to try different configurations without having to recompile. Every process will write its result to a text file in a folder called "out" so make sure you have this folder set-up. To run the code on 4 cores, use the following command: <br />
`mpiexec -n  4 Exectuable.exe` <br />
The post-processing script also needs to access to User_Input.txt so please make sure it is in the same folder or change to path within the script. The script needs to be passed the number of mpi processes that were used to run the simulation to create the gif
 file. Please pass in the number of processes used as a command line argument like (example for 4 processes used): <br />
`python post_processing.py 4`
