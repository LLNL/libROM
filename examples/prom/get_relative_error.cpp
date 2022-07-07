//               libROM MFEM Example: Relative error calculation
// 
// Compile with: ./scripts/compile.sh -m
//
// Description:  This program calculates the relative error between two solutions.
//				 The solutions are assumed to have a VisIt readable file structure.
//				 (They are generated with the "--visit" option).
//				 Relative error is calculated as:
//				 sum(abs(FOM_sol - ROM_sol) / abs(FOM_sol))
//				 In order to keep the values from blowing up, the relative error 
//				 calculation skips entries where the abs(FOM_sol) is smaller than 
//				 the system tolerance.
// 
// Usage:		 ./get_relative_error -fp Example_linear_elastic_000000/solution.000000 -rp Example_linear_elastic_rom_000000/solution.000000

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <sstream> 
#include <iterator>
#include <cassert>
#include <math.h> 
#include <limits>

using namespace std;
using namespace mfem;

// Utility functions
void load_data_into_array(double** array, string file_name);
int get_num_rows(string file_name);
double rel_error(double** fom_array, double** rom_array, int l, int d);



int main(int argc, char* argv[]) {


	// 1. Command line options
	const char *_file_name_fom = "";
	const char *_file_name_rom = "";
	int d = 2;
	int offset = 5;


	OptionsParser args(argc, argv);
	args.AddOption(&_file_name_fom, "-fp", "--fom-path", 
		"Set FOM solution path.");
	args.AddOption(&_file_name_rom, "-rp", "--rom-path", 
		"Set ROM solution path.");
	args.AddOption(&d, "-d", "--dimension",
		"Set the dimension.");
	args.AddOption(&offset, "-o", "--offset",
		"Specifies length of header in solution file. Default is 5, which works with the MFEM-VisIt file system.");
	args.Parse();

	args.PrintOptions(cout);

	// 2. Convert filenames to strings
	string file_name_fom(_file_name_fom);
	string file_name_rom(_file_name_rom);


	// 3. Get different file lengths and assert that they're equal
	int l_fom = get_num_rows(file_name_fom);
	l_fom = l_fom - offset;

	int l_rom = get_num_rows(file_name_rom);
	l_rom = l_rom - offset;

	assert(l_fom == l_rom);


	// 4. Allocate data arrays
	double** data_array_fom;
	data_array_fom = new double* [l_fom];

	double** data_array_rom;
	data_array_rom = new double* [l_rom];

	for (int i = 0; i < l_fom; i++) {
		data_array_fom[i] = new double[d];
		data_array_rom[i] = new double[d];
	}


	// 5. Load data into both arrays
	load_data_into_array(data_array_fom, file_name_fom);
	load_data_into_array(data_array_rom, file_name_rom);



	 // 6. Compute error
	double error = rel_error(data_array_fom, data_array_rom, l_fom, d);

	cout << "Error is: " << error << "\n";

	// 7. Deallocate array
	for (int i = 0; i < l_fom; i++) {
		delete[] data_array_fom[i];
		delete[] data_array_rom[i];
	}

	// 8. Delete array pointer
	delete[] data_array_fom;
	delete[] data_array_rom;

	return 1;
}



int get_num_rows(string file_name)
{
	// A function to get the number of lines in the solution matrix.
	// This is used to allocate the data arrays to be filled.


	// Initialize variables
	fstream readfile;
	int ctr = 0;

	// Loop through lines in file and increase counter
	readfile.open(file_name, ios::in); 
	if (readfile.is_open()) {   
		string tp;
		while (getline(readfile, tp)) {  
			ctr++;
		}
		readfile.close();   
	}

	return ctr;
}

void load_data_into_array(double** array, string file_name)
{
	// Function to load data into an array

	// Initialize variables
	fstream readfile;
	readfile.open(file_name, ios::in); 
	const int d = 2;
	int ctr = 0;
	int l;
	double d1, d2;


	// Loop through file and add data to array
	if (readfile.is_open()) {   
		string tp;
		while (getline(readfile, tp)) {  

			if (ctr > 4)
			{
				std::istringstream stm(tp);
				if (stm >> d1 >> d2)
				{
					array[ctr - 5][0] = d1;
					array[ctr - 5][1] = d2;
				}
			}
			ctr++;
		}
		readfile.close();   
	}

	return;

}


double rel_error(double** fom_array, double** rom_array, int l, int d)
{
	// Function that loops through the two arrays and calculates their relative
	// error.

	double error = 0;

	for (int i = 0; i < l; i++) {
		for (int j = 0; j < d; j++) {
			if (std::fabs(fom_array[i][j]) > std::numeric_limits<double>::epsilon())
			{
				error += std::fabs(fom_array[i][j] - rom_array[i][j]) / std::fabs(fom_array[i][j]);
			}
			
		}
	}

	return error;
}









