/* finite difference driver program for the wave equation.

compile with: icc -Os -o fdd fdDriver.cpp

run with: ./fdd

Sergey Voronin
*/


#include <iostream>
#include <time.h>
#include <string.h>

#include "finiteDiff.cpp"

using namespace std;


//define parameters to control which cases to display
int DISPLAY_CONSTANT = 1;
int DISPLAY_VARIABLE = 1;

double func_actual(double x, double t);
STLMatrix<double> getExactTable(long L, long T, double delx, double delt);
void printToFile(STLMatrix<double> Table, double table_x_spacing, double  table_t_spacing, char *filename);

int main(){

	cout << "Starting.." << endl;

	//problem parameters
	long L = 1;				//x limit
	long T = 1;				//t limit
	double h,k,table_x_spacing,table_t_spacing;
	long table_rows, table_cols;

	//table paremeters
	table_x_spacing = 0.1;	//table x spacing
	table_t_spacing = 0.1;	//table t spacing
	table_rows = (long)ceil(L/table_x_spacing);
	table_cols = (long)ceil(T/table_t_spacing);

	//setup tables
	STLMatrix<double> TableExact(table_rows+1,table_cols+1);
	STLMatrix<double> TableDiff(table_rows+1,table_cols+1);
	STLMatrix<double> TableN(table_rows+1,table_cols+1); 	//Table N

	//setup fdDiff object
	finiteDiff fdObject;

	//setup h and k value arrays to use
	int num_test_vals = 4;
	//double h_array[4] = { .1, .025, .00625, .0015625  };
	//double k_array[4] = { .02, .005, .00125, 3.125e-4 };
	double h_array[4] = { .1, .1/4.0, .1/16.0, .1/64.0  };
	double k_array[4] = { .02, .02/4.0, .02/16.0, .02/64.0 };
	double max_error[4];

	//setup error storage containers
	STLMatrix<double> ConstantSTDErrors(num_test_vals,1);
	STLMatrix<double> ConstantNONSTDErrors(num_test_vals,1);
	STLMatrix<double> VariableSTDErrors(num_test_vals,1);
	STLMatrix<double> VariableNONSTDErrors(num_test_vals,1);



	if(DISPLAY_CONSTANT){
		/* *** CONSTANT COEFFICIENT CASE - STANDARD ALGORITHM *** */
	
		cout << "\n*** CONSTANT COEFFICIENT CASE - STANDARD ALGORITHM ***\n" << endl;
	
		//get exact table
		TableExact = getExactTable(L, T, table_x_spacing, table_t_spacing);
	
		//print the exact table to file
		cout << "writing exact table to file.." << endl;
		printToFile(TableExact, table_x_spacing, table_t_spacing, "out_constant_exact.txt");
		

		//loop through the h,k arrays and compare approximation tables to the exact one
		for(int ind=0; ind<num_test_vals; ind++){
			h = h_array[ind];
			k = k_array[ind];
	
			//initialize object w/ above parameters
			fdObject.init(L, T, h, k, 1);
	
			//grab table
			TableN = fdObject.getSTDTable(table_x_spacing, table_t_spacing);
	
			//print table
			//cout << "Table " << ind << ":\n" << TableN << endl;
	
			//get max absolute error
			TableDiff = (TableExact - TableN).abs();
			max_error[ind] = TableDiff.max();
	
			if(ind == 3){
				//print the error table to file
				printToFile(TableDiff, table_x_spacing, table_t_spacing, "out_constant_std.txt");
				
			}
	
			//record error
			ConstantSTDErrors.set(ind,0, max_error[ind]);
	
			//cout error
			cout << "h = " << h << " & k = " << k << endl;
			cout << "\n\tmax error # " << ind << " = " << max_error[ind] << "\n\n";
		}
	
	
		//calculate ratios
		cout << "Error Ratios: " << endl;
		for(int ind=1; ind<num_test_vals; ind++)
			cout << (double)(max_error[ind]/max_error[ind-1]) << "\t";
		cout << endl;
	
	
	
		/* *** CONSTANT COEFFICIENT CASE - NON-STANDARD ALGORITHM *** */
	
		cout << "\n\n*** CONSTANT COEFFICIENT CASE - NON-STANDARD ALGORITHM ***\n" << endl;
	
		//get exact table
		TableExact = getExactTable(L, T, table_x_spacing, table_t_spacing);

	
		//loop through the h,k arrays and compare approximation tables to the exact one
		for(int ind=0; ind<num_test_vals; ind++){
			h = h_array[ind];
			k = k_array[ind];
	
			//initialize object w/ above parameters
			fdObject.init(L, T, h, k, 1);
	
			//grab table
			TableN = fdObject.getNONSTDTable(table_x_spacing, table_t_spacing);
	
			//print table
			//cout << "Table " << ind << ":\n" << TableN << endl;
	
			//get max absolute error
			TableDiff = (TableExact - TableN).abs();
			max_error[ind] = TableDiff.max();
	

			if(ind == 3){
				//print the error table to file
				printToFile(TableDiff, table_x_spacing, table_t_spacing, "out_constant_nonstd.txt");
				
			}

			
			//record error
			ConstantNONSTDErrors.set(ind,0, max_error[ind]);
	
			//cout error
			cout << "h = " << h << " & k = " << k << endl;
			cout << "\n\tmax error # " << ind << " = " << max_error[ind] << "\n\n";
		}
	
		//calculate ratios
		cout << "Error Ratios: " << endl;
		for(int ind=1; ind<num_test_vals; ind++)
			cout << (double)(max_error[ind]/max_error[ind-1]) << "\t";
		cout << endl;

	}



	if(DISPLAY_VARIABLE){
		/* *** VARIABLE COEFFICIENT CASE - STANDARD ALGORITHM *** */
	
		cout << "\n\n*** VARIABLE COEFFICIENT CASE - STANDARD ALGORITHM ***\n" << endl;
	
		//get "exact" table - in this case, we must use a fine approximation
		cout << "setting up exact table.. (this may take a while)" << endl;
		h = 0.001/2.0;
		k = 0.0002/2.0;
		fdObject.init(L, T, h, k, 0);
		TableExact = fdObject.getSTDTable(table_x_spacing, table_t_spacing);
	
		//print the exact table to file
		cout << "writing exact table to file.." << endl;
		printToFile(TableExact, table_x_spacing, table_t_spacing, "out_variable_exact.txt");
		

		//loop through the h,k arrays and compare approximation tables to the exact one
		for(int ind=0; ind<num_test_vals; ind++){
			h = h_array[ind];
			k = k_array[ind];
	
			//initialize object w/ above parameters
			fdObject.init(L, T, h, k, 0);
	
			//grab table
			TableN = fdObject.getSTDTable(table_x_spacing, table_t_spacing);
	
			//print table
			//cout << "Table " << ind << ":\n" << TableN << endl;
	
			//get max absolute error
			TableDiff = (TableExact - TableN).abs();
			max_error[ind] = TableDiff.max();


			if(ind == 3){
				//print the error table to file
				printToFile(TableDiff, table_x_spacing, table_t_spacing, "out_variable_std.txt");
				
			}
			
			
			//record error
			VariableSTDErrors.set(ind,0, max_error[ind]);
	
			//cout error
			cout << "h = " << h << " & k = " << k << endl;
			cout << "\n\tmax error # " << ind << " = " << max_error[ind] << "\n\n";
		}
	
	
		//calculate ratios
		cout << "Error Ratios: " << endl;
		for(int ind=1; ind<num_test_vals; ind++)
			cout << (double)(max_error[ind]/max_error[ind-1]) << "\t";
		cout << endl;
	
	
	
		/* *** VARIABLE COEFFICIENT CASE - NON-STANDARD ALGORITHM *** */
	
		cout << "\n\n*** VARIABLE COEFFICIENT CASE - NON-STANDARD ALGORITHM ***\n" << endl;
	
		//get "exact" table - in this case, we must use a fine approximation
		cout << "setting up exact (fine approx) table.. (this may take a while)" << endl;
		h = 0.001/2.0;
		k = 0.0002/2.0;
		fdObject.init(L, T, h, k, 0);
		TableExact = fdObject.getNONSTDTable(table_x_spacing, table_t_spacing);
	
	
		//loop through the h,k arrays and compare approximation tables to the exact one
		for(int ind=0; ind<num_test_vals; ind++){
			h = h_array[ind];
			k = k_array[ind];
	
			//initialize object w/ above parameters
			fdObject.init(L, T, h, k, 0);
	
			//grab table
			TableN = fdObject.getNONSTDTable(table_x_spacing, table_t_spacing);
	
			//print table
			//cout << "Table " << ind << ":\n" << TableN << endl;
	
			//get max absolute error
			TableDiff = (TableExact - TableN).abs();
			max_error[ind] = TableDiff.max();
	

			if(ind == 3){
				//print the error table to file
				printToFile(TableDiff, table_x_spacing, table_t_spacing, "out_variable_nonstd.txt");
				
			}
			
			
			//record error
			VariableNONSTDErrors.set(ind,0, max_error[ind]);
	
			//cout error
			cout << "h = " << h << " & k = " << k << endl;
			cout << "\n\tmax error # " << ind << " = " << max_error[ind] << "\n\n";
		}
	
	
		//calculate ratios
		cout << "Error Ratios: " << endl;
		for(int ind=1; ind<num_test_vals; ind++)
			cout << (double)(max_error[ind]/max_error[ind-1]) << "\t";
		cout << endl;
	
		//end variable
	}


	//print summary of errors
	cout << "\n\nSummary of Errors for Our Problem:\n" << endl;

	cout << "f(x) = sin(" << M << "*pi*x)" << "   g(x) = 0" << "   v_constant = " << V << "   v_variable(x) = 1 + x^2" << "   w = " << M << "*pi" << endl;
	
	if(DISPLAY_CONSTANT){
		cout << "\n*** Constant Coefficient Errors (STD & NON-STD) ***\n" << endl;
		//cout << ConstantSTDErrors << endl;
		//cout << ConstantNONSTDErrors << endl;
		cout << "Standard Algorithm Errors: " << endl;
		for(int ind=0; ind<num_test_vals; ind++){
			cout << "h = " << h_array[ind] << "\tk = " << k_array[ind] << "\tmax absolute error = " << ConstantSTDErrors.get(ind,0) << endl;
		}

		cout << "\nNon-Standard Algorithm Errors: " << endl;
		for(int ind=0; ind<num_test_vals; ind++){
			cout << "h = " << h_array[ind] << "\tk = " << k_array[ind] << "\tmax absolute error = " << ConstantNONSTDErrors.get(ind,0) << endl;
		}
	}

	if(DISPLAY_VARIABLE){
		cout << "\n\n*** Variable Coefficient Errors (STD & NON-STD) ***\n" << endl;
		//cout << VariableSTDErrors << endl;
		//cout << VariableNONSTDErrors << endl;
		//cout << "\n\n";

		cout << "Standard Algorithm Errors: " << endl;
		for(int ind=0; ind<num_test_vals; ind++){
			cout << "h = " << h_array[ind] << "\tk = " << k_array[ind] << "\tmax absolute error = " << VariableSTDErrors.get(ind,0) << endl;
		}

		cout << "\nNon-Standard Algorithm Errors: " << endl;
		for(int ind=0; ind<num_test_vals; ind++){
			cout << "h = " << h_array[ind] << "\tk = " << k_array[ind] << "\tmax absolute error = " << VariableNONSTDErrors.get(ind,0) << endl;
		}
	}


	//print processing time
	time_t time = clock();
	cout << "\n\nprocessor time used: " << (double)(time*1.0/CLOCKS_PER_SEC) << " sec " << endl;


	//return and exit
	return 0;
}



/* exact solution for comparison in constant coefficient case */
double func_actual(double x, double t){
	double res;
	double pi = M_PI;
	double m = M;
	double v = V;
	res = sin(m*pi*x)*cos(v*m*pi*t);
	return res;
}



/* returns a table with the exact solution (for the constant coefficient case) */
STLMatrix<double> getExactTable(long L, long T, double delx, double delt){

	long i,j;
	double x,t,val;
	double h = delx;
	double k = delt;
	long num_i = (long)ceil(L/h);
	long num_j = (long)ceil(T/k);
	STLMatrix<double> Table(num_i+1,num_j+1);

	for(j=0; j<=num_j; j++){
		for(i=0; i<=num_i; i++){
			if(i==0 || i==num_i)
				val = 0;
			else{
				x = i*h;
				t = j*k;
				val = func_actual(x,t);
			}
			Table.set(i,j, val);
		}
	}
	
	return Table.transpose();
}


/* prints a table to a file */
void printToFile(STLMatrix<double> Table, double table_x_spacing, double table_t_spacing, char *filename){

	char abs_filename[30];
	strcpy(abs_filename,"data_files/");
	strcat(abs_filename,filename);

	fstream ftable(abs_filename,ios::out);
	long num_rows = Table.getNumRows();
	long num_cols = Table.getNumCols();
	long row,col;
	double x,t,val;

	for(col=0; col<num_cols; col++){
		for(row=0; row<num_rows; row++){
			val = Table.get(row,col);
			x = table_x_spacing*row;
			t = table_t_spacing*col;
			ftable << x << "\t" << t << "\t" << val << endl;
		}
	}

	ftable.close();
}


