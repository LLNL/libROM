// Tutorial from https://www.tutorialspoint.com/read-data-from-a-text-file-using-cplusplus
// double reading from https://stackoverflow.com/questions/25444449/how-do-i-convert-a-stdstring-containing-doubles-to-a-vector-of-doubles
//string file_name_fom = "Example_linear_elastic_000000/solution.000000";
//string file_name_rom = "Example_linear_elastic_rom_000000/solution.000000";


#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <sstream> 
#include <iterator>
#include <cassert>
#include <math.h> 

using namespace std;
using namespace mfem;

void load_data_into_array(double** array, string file_name);
int get_num_rows(string file_name);
double rel_error(double** fom_array, double** rom_array, int l, int d);



int main(int argc, char* argv[]) {




	

	string file_name_fom;
	string file_name_rom;

	int d = 2;
	int offset = 5;


	OptionsParser args(argc, argv);
	args.AddOption(&file_name_fom, "-fp", "--fom-path", 
		"Set FOM solution path.");
	args.AddOption(&file_name_rom, "-rp", "--rom-path", 
		"Set ROM solution path.");
	args.AddOption(&d, "-d", "--dimension",
		"Set the dimension.");
	args.Parse();

	// Get different file lengths and assert that they're equal
	int l_fom = get_num_rows(file_name_fom);
	l_fom = l_fom - offset;

	int l_rom = get_num_rows(file_name_rom);
	l_rom = l_rom - offset;

	assert(l_fom == l_rom);


	// Allocate data array1 and 2
	double** data_array_fom;
	data_array_fom = new double* [l_fom];

	double** data_array_rom;
	data_array_rom = new double* [l_rom];

	for (int i = 0; i < l_fom; i++) {
		data_array_fom[i] = new double[d];
		data_array_rom[i] = new double[d];
	}


	// Load data into both arrays
	load_data_into_array(data_array_fom, file_name_fom);
	load_data_into_array(data_array_rom, file_name_rom);


	/*
	for (int i = 0; i < l_fom; i++) {

		cout << data_array_fom[i][0] << " " << data_array_fom[i][0] << "\n";

	}
	 */

	 // Compute error
	double error = rel_error(data_array_fom, data_array_rom, l_fom, d);

	cout << "Error is: " << error << "\n";

	// Deallocate array
	for (int i = 0; i < l_fom; i++) {
		delete[] data_array_fom[i];
		delete[] data_array_rom[i];
	}

	// Delete array pointer
	delete[] data_array_fom;
	delete[] data_array_rom;

	return 1;
}



int get_num_rows(string file_name)
{
	fstream newfile;
	int ctr = 0;
	newfile.open(file_name, ios::in); //open a file to perform read operation using file object
	if (newfile.is_open()) {   //checking whether the file is open
		string tp;
		while (getline(newfile, tp)) {  //read data from file object and put it into string.
			ctr++;
		}
		newfile.close();   //close the file object.
	}

	return ctr;
}

void load_data_into_array(double** array, string file_name)
{
	fstream newfile;
	newfile.open(file_name, ios::in); //open a file to perform read operation using file object
	const int d = 2;
	int ctr = 0;
	int l;
	double d1, d2;

	if (newfile.is_open()) {   //checking whether the file is open
	//cout << "All good!\n";
		string tp;
		while (getline(newfile, tp)) {  //read data from file object and put it into string.

			if (ctr > 4)
			{
				// create an input string stream to read from the string
				std::istringstream stm(tp);
				if (stm >> d1 >> d2)
				{
					array[ctr - 5][0] = d1;
					array[ctr - 5][1] = d2;
				}


			}

			ctr++;
		}

		//cout << "Length of list is: " << l;
		newfile.close();   //close the file object.
	}

	return;

}


double rel_error(double** fom_array, double** rom_array, int l, int d)
{
	double error = 0;
	//double tol = 1.0e0;

	for (int i = 0; i < l; i++) {
		for (int j = 0; j < d; j++) {

			error += std::fabs(fom_array[i][j] - rom_array[i][j]) /// (std::fabs(fom_array[i][j]) + tol);
		}
	}

	return error;
}









