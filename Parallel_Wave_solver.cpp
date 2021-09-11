#define _USE_MATH_DEFINES

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <chrono>
#include <omp.h>


using namespace std;


int id, p;  // process id and global number of processes  used
int* start_row; // array of start row number for every process
int neigbour_up; // store id of the neigbour above 
int neigbour_down; // store id of the neigbour below 
int cnt; // count number of mpi operations 
double y_proc_start; // y value from 0 where the process first row is located
int remaining_rows; // rows left during set-up 

int periodic_BC, non_periodic_BC; // store choice of B.C. 

double x_disturbance;
double y_disturbance;
double y_barrier;
double r_distrubance;


MPI_Request* request; // MPI request array to store mpi operations

double* recv_buffer_up; // store values send from process above
double* recv_buffer_down; // store values send from process below


// vector grid was replaced with array which has improved the speed 
double* grid;
double* new_grid;
double* old_grid;
int imax, jmax; //numbers of grid points in y and x
double t_max; //stopping t value
double t, t_out = 0.0, dt_out, dt; // time, timings for output, timesteps for output, simulation timestep
double y_max, x_max, dx, dy; // domain y, domain x, dx, dx
double c; // speed of wave

// dataypes for two different mpi operations
MPI_Datatype Datatype_top, Datatype_bottom;

void createdatatypes(double* data, int imax, int jmax)
{	// Function to set up the MPI datatypes for the send
	// Send top row and bottom row of array

	vector<int> block_lengths; // length of array in y
	vector<MPI_Datatype> typelist; // type 
	vector<MPI_Aint> addresses; // memory address of array
	MPI_Aint add_start; // memory address of first point that needs to be send


	MPI_Get_address(data, &add_start);


	//Send top row of array
	int block_length = jmax;
	MPI_Datatype typeval = MPI_DOUBLE;
	MPI_Aint address;
	MPI_Get_address(&data[0], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_top);
	MPI_Type_commit(&Datatype_top);

	//send bottom row of array
	MPI_Get_address(&data[jmax*(imax - 1)], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_bottom);
	MPI_Type_commit(&Datatype_bottom);
}


void grid_to_file(int out)
{	//Store the current itteration in a txt file

	
	stringstream fname;
	fstream f1;
	//name of the output file uses if and output number
	//this is important for post-processing  script 
	fname << "./out/output" << "_" << id << "_" << out << ".txt";
	f1.open(fname.str().c_str(), ios_base::out);

	//loop over all rows
	for (int i = 0; i < imax; i++)
	{	//loop over all columns 
		for (int j = 0; j < jmax; j++)
			f1 << grid[i * jmax + j] << "\t";
		f1 << endl;
	}
	f1.close();
}




int read_send_inputfile(void)
{	  //The process with ID=0 reads some user defined varibles 
	  //from User_Input.txt and sends it to all other processes 
	  //via a MPI blocking communication. All processes will store the 
	  //recived values in the right varibles	
	
	int file_not_exit = 0;

	double *save_input = new double[13]; //array to store data from file
	if (id == 0)
	{
		
		//open User_Input.txt file
		const string name_file = "User_Input.txt";
		ifstream file(name_file);
		// dumy varible to store data from file
		string line;

		//if file is not in the folder, print out and return a 1 to main
		//so that the file can be stoped from main
		if (file.good()) {

			// count number of lines for correct postion in the array
			int line_cnt = 0;
			while (!file.eof()) {
				//get every line
				getline(file, line);

				// if a line does not start if / we need to save it
				if (line[0] != '/') {
					//store every line in array as double
					save_input[line_cnt] = stod(line.c_str());
					line_cnt++;
				}

			}


		}
		else
		{
			// In case the file does not exit, a Bcast needs to happen for other processes to continue otherwise it will get stuck
			MPI_Bcast(save_input, 13, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// print message to screen and exit with error code 1 to main
			cout << "No User_Input.txt in the folder, please copy file from Github \n";
			return 1;
		

		}


	
	}

	
	//send data to other processes
	MPI_Bcast(save_input, 13, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	

	//convert the values stored in the array if necessary and 
	//set it to the correct variable 
	y_max = save_input[0];
	x_max = save_input[1];
	imax = (int)save_input[2];
	jmax = (int)save_input[3];
	t_max = save_input[4];
	dt_out = save_input[5];
	c = save_input[6];
	r_distrubance = save_input[7];
	x_disturbance = save_input[8];
	y_disturbance = save_input[9];
	y_barrier = save_input[10];
	periodic_BC = (int)save_input[11];
	non_periodic_BC = (int)save_input[12];

	//if everything has finished return a 0 so that main knows that all
	//variable are set and simulation can start 
	return 0;

}



void inner_itteration(void)
{	//This does the inner itteration but does not include boundary row on the left/right and top/bottom

	for (int i = 1; i < imax-1; i++)
		for (int j = 1; j < jmax-1; j++)
			new_grid[i * jmax + j] = pow(dt * c, 2.0) * ((grid[(i+1) * jmax + j] - 2.0 * grid[i * jmax + j] + grid[(i - 1) * jmax + j]) / pow(dy, 2.0) + (grid[i * jmax + j + 1] - 2.0 * grid[i * jmax + j] + grid[i * jmax + j - 1]) / pow(dx, 2.0)) + 2.0 * grid[i * jmax + j] - old_grid[i * jmax + j];


}


void do_itteration()
{
	// set mpi counter to zero
	cnt = 0;

	

	// Do communication between processes, first set up recv and then send
	// Store the recv in a buffer array which has the same size
	MPI_Irecv(&recv_buffer_up[0], jmax, MPI_DOUBLE, neigbour_up, 0, MPI_COMM_WORLD, &request[cnt]);
	cnt++;
	MPI_Irecv(&recv_buffer_down[0], jmax , MPI_DOUBLE, neigbour_down, 1, MPI_COMM_WORLD, &request[cnt]);
	cnt++;
	MPI_Isend(grid, 1, Datatype_top, neigbour_up, 1, MPI_COMM_WORLD, &request[cnt]);
	cnt++;
	MPI_Isend(grid, 1, Datatype_bottom, neigbour_down, 0, MPI_COMM_WORLD, &request[cnt]);
	cnt++;
	
	

	// while the mpi comms is working the the background, the itteration in the middle of the domain is done
	inner_itteration();


	// Once the itteration in the domain is finished, we need to wait for all comms to finish such that every proccess
	// has the infomation it needs to do the boundary of its domain
	MPI_Waitall(cnt, request, MPI_STATUSES_IGNORE);

	// The B.C. depend on user input
	if (periodic_BC == 1)
	{
		
		// do the periodic B.C on the right and left hand side of the domain by using points from the opposite side of domain
		for (int i = 1; i < imax - 1; i++)
		{

			new_grid[jmax*i] = pow(dt * c, 2.0) * ((grid[jmax*(i+1)] - 2.0 * grid[jmax*i] + grid[jmax*(i-1)]) / pow(dy, 2.0) +
				(grid[jmax*i+1] - 2.0 * grid[jmax*i] + grid[jmax*i+(jmax-1)]) / pow(dx, 2.0)) + 2.0 * grid[jmax*i] - old_grid[jmax*i];

			new_grid[jmax*i+(jmax-1)] = pow(dt * c, 2.0) * ((grid[jmax*(i+1)+jmax-1] - 2.0 * grid[jmax*i+jmax-1] + grid[jmax*(i-1)+jmax-1]) / pow(dy, 2.0) 
				+ (grid[jmax*i] - 2.0 * grid[jmax*i+jmax-1] + grid[jmax*i+jmax-2]) / pow(dx, 2.0)) + 2.0 * grid[jmax*i+jmax-1] - old_grid[jmax*i+jmax-1];
		}
		
		
		
		// This implements the top and bottom rows which uses points from other processes via the recv_buffer arrays
		for (int j = 1; j < jmax - 1; j++)
		{
			// top row
			new_grid[j] = pow(dt * c, 2.0) * ((grid[jmax + j] - 2.0 * grid[j] + recv_buffer_up[j]) / pow(dy, 2.0) + 
				(grid[j+1] - 2.0 * grid[j] + grid[j-1]) / pow(dx, 2.0)) + 2.0 * grid[j] - old_grid[j];

			// bottom row
			new_grid[jmax*(imax-1)+j] = pow(dt * c, 2.0) * ((recv_buffer_down[j] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-2)+j]) / pow(dy, 2.0) + 
				(grid[jmax*(imax-1)+j+1] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-1)+j-1]) / pow(dx, 2.0)) + 2.0 * grid[jmax*(imax-1)+j] - old_grid[jmax*(imax-1)+j];
		}

		// Also do this for top left
		new_grid[0] = pow(dt * c, 2.0) * ((grid[jmax] - 2.0 * grid[0] + recv_buffer_up[0]) / pow(dy, 2.0) +
			(grid[1] - 2.0 * grid[0] + grid[jmax - 1]) / pow(dx, 2.0)) + 2.0 * grid[0] - old_grid[0];
		// and bottom left corner
		new_grid[jmax * (imax - 1)] = pow(dt * c, 2.0) * ((recv_buffer_down[0] - 2.0 * grid[jmax * (imax - 1)] + grid[jmax * (imax - 2)]) / pow(dy, 2.0) +
			(grid[jmax * (imax - 1) + 1] - 2.0 * grid[jmax * (imax - 1)] + grid[jmax * (imax - 1) + jmax - 1]) / pow(dx, 2.0)) + 2.0 * grid[jmax * (imax - 1)] - old_grid[jmax * (imax - 1)];

		// and for top right and 
		new_grid[jmax-1] = pow(dt * c, 2.0) * ((grid[jmax + jmax - 1] - 2.0 * grid[jmax-1] + recv_buffer_up[jmax - 1]) / pow(dy, 2.0) + 
			(grid[0] - 2.0 * grid[jmax-1] + grid[jmax - 2]) / pow(dx, 2.0)) + 2.0 * grid[jmax-1] - old_grid[jmax-1];
		// bottom right corner
		new_grid[jmax*(imax-1)+jmax-1] = pow(dt * c, 2.0) * ((recv_buffer_down[jmax - 1] - 2.0 * grid[jmax*(imax-1)+jmax-1] + grid[jmax*(imax-2)+jmax-1]) / pow(dy, 2.0) +
			(grid[jmax*(imax-1)] - 2.0 * grid[jmax*(imax-1)+jmax-1] + grid[jmax*(imax-1)+jmax-2]) / pow(dx, 2.0)) + 2.0 * grid[jmax*(imax-1)+jmax-1] - old_grid[jmax*(imax-1)+jmax-1];
		
	}
	// This are the non peridoic B.C
	else
	{	// Neumann 
		if (non_periodic_BC == 0)
		{

			// do top and bottom row which uses points from other processes via the recv_buffer arrays
			for (int j = 1; j < jmax - 1; j++)
			{
				new_grid[j] = pow(dt * c, 2.0) * ((grid[jmax + j] - 2.0 * grid[j] + recv_buffer_up[j]) / pow(dy, 2.0) + (grid[j+1] - 2.0 * grid[j] + grid[j-1]) / pow(dx, 2.0)) + 2.0 * grid[j] - old_grid[j];
				new_grid[jmax*(imax-1)+j] = pow(dt * c, 2.0) * ((recv_buffer_down[j] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-2)+j]) / pow(dy, 2.0) + (grid[jmax*(imax-1)+j+1] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-1)+j-1]) / pow(dx, 2.0)) + 2.0 * grid[jmax*(imax-1)+j] - old_grid[jmax*(imax-1)+j];
			}

			// left and right hand corner of array
			for (int i = 0; i < imax; i++)
			{

				new_grid[jmax * i] = new_grid[jmax * i + 1];
				new_grid[jmax * i + jmax - 1] = new_grid[jmax * i + jmax - 2];
			}


			// If the process if at the top or at the bottom of the global domain i.e. id = 1 or id = p - 1,
			// it needs to enfornce the B.C. instead of the using the values above
			if (id == 0)
			{
				// set bottom row 
				for (int j = 0; j < jmax; j++)
				{
					new_grid[jmax*(imax-1)+j] = new_grid[jmax*(imax-2)+j];
				}
				

			}

			if (id == p - 1)
			{
				// set top row
				for (int j = 0; j < jmax; j++)
				{
					new_grid[j] = new_grid[jmax + j];
					
				}


			}
		


		}
		// dirichlet B.C.
		else {


			// do top and bottom row which uses points from other processes via the recv_buffer arrays
			for (int j = 1; j < jmax - 1; j++)
			{
				new_grid[j] = pow(dt * c, 2.0) * ((grid[jmax + j] - 2.0 * grid[j] + recv_buffer_up[j]) / pow(dy, 2.0) + (grid[j+1] - 2.0 * grid[j] + grid[j-1]) / pow(dx, 2.0)) + 2.0 * grid[j] - old_grid[j];
				new_grid[jmax*(imax-1)+j] = pow(dt * c, 2.0) * ((recv_buffer_down[j] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-2)+j]) / pow(dy, 2.0) + (grid[jmax*(imax-1)+j+1] - 2.0 * grid[jmax*(imax-1)+j] + grid[jmax*(imax-1)+j-1]) / pow(dx, 2.0)) + 2.0 * grid[jmax*(imax-1)+j] - old_grid[jmax*(imax-1)+j];
			}
					

			// left and right hand corner of array
			for (int i = 0; i < imax; i++)
			{

				new_grid[jmax * i] = 0;
				new_grid[jmax * i + jmax - 1] = 0;
			}


			// If the process if at the top or at the bottom of the global domain i.e. id = 1 or id = p - 1,
			// it needs to enfornce the B.C. instead of the using the values above
			if (id == 0)
			{
				for (int j = 0; j < jmax; j++)
				{
					new_grid[jmax*(imax-1)+j] = 0;
				}


			}

			if (id == p - 1)
			{
				for (int j = 0; j < jmax; j++)
				{
					new_grid[j] = 0;
				}


			}
			

		}
	}


	// if y_barrier is set in User_input file
	if (y_barrier != 0) {


		//check if start and end y value of process is between y_barrier
		if (start_row[id] * dy < y_barrier && start_row[id + 1] * dy > y_barrier) {


			// find row where barrier is set
			int row = (int)((imax - (y_barrier/dy - start_row[id]))*dy / dy);


			// set barrier from 0 to half the domain
			for (int j = 0; j < (jmax / 2); j++)
			{
				new_grid[jmax * row + j] = 0;

			}
		}
	}
	
	
	// advance time by dt
	t += dt;

	// swap pointer such that old_grid = grid, grid = ned_grid and new_grid = old_grid
	swap(old_grid, new_grid); 
	swap(old_grid, grid); 



	
}


void set_up_parameters() {
	// This function sets up the rectangular domain for every process
	// Every process has the same number of columns but different number of rows

	// dy and dx are the same in whole domain so can be done using the y_max/x_max values
	dx = x_max / ((double)jmax - 1);
	dy = y_max / ((double)imax - 1);

	// start_row is array to store the first for every process. For exmaple imax = 100 and 4 processes:
	// start_row = [0, 25, 50, 75, 100, 100]
	// start_row[id] is the starting row for a process
	// start_row[id + 1] - 1 is the end row for a process
	// start_row[p + 1] is the overall number of rows
	start_row = new int[p + 2];
	start_row[0] = 0;
	start_row[p+1] = imax;
	start_row[p] = imax;

	remaining_rows = imax;

	// this is only done if more than one process is used
	if (p != 1) {

		// Find start_row for every process
		for (int i = 0; i < p; i++)
		{
			int nr_rows_process = remaining_rows / (p - i);

			if (i == 0)
			{
				start_row[i + 1] = nr_rows_process;
			}


			if (i + 1 < p && i != 0)
			{
				start_row[i + 1] = nr_rows_process + start_row[i];

			}

			remaining_rows -= nr_rows_process;
		}

		// update imax for every process, new imax depends on number of rows for every process
		if (id != p)
		{
			imax = start_row[id + 1] - start_row[id];

		}
		else
		{
			imax = start_row[p] - start_row[id - 1];

		}

		// store id of neigbours
		neigbour_up = id + 1;
		neigbour_down = id - 1;


		// if id = 0, neigbour below is process on top of array
		if (id == 0)
		{
			neigbour_up = id + 1;
			neigbour_down = p - 1;

		}

		// if id = p - 1, neigbour below is process on bottom of array
		if (id == p - 1)
		{
			neigbour_up = 0;
			neigbour_down = p - 2;

		}

		// y_proc_start is the y value of the first row of the process inside the global array
		if (id != p)
		{
			y_proc_start = dy * start_row[id];

		}
		else
		{
			y_proc_start = dy * start_row[id - 1];
		}

	}

	// since we now know the size of our local arrays, we can set up the grid and buffer 
	request = new MPI_Request[2 * 2];
	recv_buffer_up = new double[jmax];
	recv_buffer_down = new double[jmax];
	grid = new double[imax * jmax];
	new_grid = new double[imax * jmax];
	old_grid = new double[imax * jmax];

}


void inital_condition() {
	// This addes a inital condition to the array


	double r_splash = r_distrubance; // set radius of disturbance 
	double x_splash = x_disturbance; // set x posision
	double y_splash = y_disturbance; // set y position
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			// determine global x and y values for every point in the local arrays
			double x = dx * j;
			double y = (y_proc_start + dy * (start_row[id + 1] - start_row[id])) - dy * i;

			double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

			// set the grid to zero to override values that might be stored in there
			grid[i * jmax + j] = 0;
			old_grid[i * jmax + j] = 0;

			
			if (dist < r_splash)
			{
				double h = 5.0 * (cos(dist / r_splash * M_PI) + 1.0);

				grid[i * jmax + j] = h;
				old_grid[i * jmax + j] = h;
			}
		}
}

int main(int argc, char *argv[]) {
	
	// start time
	auto start_time = chrono::high_resolution_clock::now();

	// MPI values
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	
	int status_; // stores value of read_send_inputfile function
	status_ = read_send_inputfile(); // read user inout via User_Input.txt

	// if read_send_inputfile is not successful it will return a 1 which means that 
	// the main needs to stopped by calling MPI_Finalize
	if (status_ == 1)
	{
		cout << "Stopping, please restart executable with User_Input file in the same folder";
		MPI_Finalize();
		return 1;
	}

	
	// set up local parameters
	set_up_parameters();
	// set up MPI data types
	createdatatypes(grid, imax, jmax);


	

	t = 0.0;
	dt = 0.4 * min(dx, dy) / c;

	int out_cnt = 0;
	int it = 0;
	
	// set up intial disturbance on grid
	inital_condition();


	// save intial condion to text file
	grid_to_file(out_cnt);
	out_cnt++;
	t_out += dt_out;


	

	
	// loop until t_max is reached
	while (t < t_max)
	{
		// do itteration 
		do_itteration();
		
		// if it is time to save grid to file
		if (t_out <= t)
		{
			// Only one process should be print out an update to see if the code is running as expected
			if (id == 0)
			{ 
				cout << "output: " << out_cnt << "\tt: " << t << "\titeration: " << it << endl;
			}
			// save grid to file
			grid_to_file(out_cnt);
			out_cnt++;
			t_out += dt_out;
		}

		it++;
	}
	

	// once everything it done call MPI_Finalize
	MPI_Finalize();

	// delete arrays
	delete[] request;
	delete[] recv_buffer_down;
	delete[] recv_buffer_up;
	
	// Timings for finish
	auto end_time = chrono::high_resolution_clock::now(); //end time
	chrono::duration<double> elapsed = end_time - start_time;

	// process 0 prints out the time it has taken
	if (id == 0)
	{
		std::cout << "\nTime: " << elapsed.count(); // print total time taken

	}


	return 0;
}
