/* finite difference class for the wave equation 

Sergey Voronin
*/

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <time.h>

//load STLMatrix
#include "STLMatrix.cpp"

//define globar input angular frequency number ( f = sin(M*pi*x) )
double M = 8.0;
//define global constant velocity for constant coeff problem
double V = 2.0;


using namespace std;

class finiteDiff {

	private:

		//general parameters --------------------------->

		//initial conditions and settings
		long L;					//x limit
		long T;					//t limit
		double h;				//x spacing
		double k;				//t spacing
		double tol;				//numbers less than tol set to 0
		int constant_coeff; 	//whether or not this is a constant coeff problem
		
		//table variables
		double table_x_spacing;
		double table_t_spacing;



		//class variables ------------------------->
		double val,actual_val,abs_error,x,t,temp;
		double con;
		long i,j;

		//functions ------------------------------->
		double func_f(double x);
		double func_g(double x);
		double func_v(double x);
		double func_v_constant(double x);
		double func_u(double x);
		double func_u_constant(double x);


	public:
	
		/* constructor */
		finiteDiff();
		
		/* utility functions */
		void init(long L_val, long T_val, double h_val, double k_val, int constant_coeff_val);
		STLMatrix<double> getSTDTable(double table_x_spacing_val, double table_t_spacing_val);
		STLMatrix<double> getNONSTDTable(double table_x_spacing_val, double table_t_spacing_val);
};


/* default constructor*/
finiteDiff::finiteDiff(){}


/* reinitializes object using the new parameters */
void finiteDiff::init(long L_val, long T_val, double h_val, double k_val, int constant_coeff_val){
	L = L_val;
	T = T_val;
	h = h_val;
	k = k_val;
	constant_coeff = constant_coeff_val;

	tol = 1e-12;
	con = (k*k)/(h*h);
}


/* returns a table of approximated values via finite differences (std algorithm) */
STLMatrix<double> finiteDiff::getSTDTable(double table_x_spacing_val, double table_t_spacing_val){

	table_x_spacing = table_x_spacing_val;
	table_t_spacing = table_t_spacing_val;

	//setup finite difference variables and storage containers
	long num_i = (long)ceil(L/h);
	long num_j = (long)ceil(T/k);

	cout << "setting up structures.. " << endl;

	//setup storage containers for the 3 rows we need
	STLMatrix<double> rowj(num_i+1,1);		//row j;
	STLMatrix<double> rowj1(num_i+1,1);		//row j-1;
	STLMatrix<double> row_temp(num_i+1,1);	//temp row

	//setup table variables
	long delx = (long)ceil(table_x_spacing/h);
	long delt = (long)ceil(table_t_spacing/k);
	long table_rows = (long)ceil(L/table_x_spacing);
	long table_cols = (long)ceil(T/table_t_spacing);
	long table_row = 0, table_col = 0;
	int recorded_entry = 0;
	STLMatrix<double> Table(table_rows+1,table_cols+1);
	

	cout << "starting finite differences.. " << endl;

	//perform finite differences -------------------->
	
	//step1 -> U(i,0) = f(ih) so fill up row(j-1)
	cout << "U(i,0) ------> \n\n";
	recorded_entry = 0;
	j = 0;
	for(i=0; i<=num_i; i++){
		if(i==0 || i==num_i)
			val = 0;
		else{
			val = func_f(i*h);
			if(abs(val)<tol)
				val = 0;
		}
		
		//set row j-1 to val at this i
		rowj1.set(i,0, val);

		//possibly print table entry
		if( fmod((double)i,(double)delx)<tol && fmod((double)j,(double)delt)<tol ){
			Table.set(table_row,table_col, val);
			table_col++;
			recorded_entry = 1;
		}
	}
	//update table row and col
	if(recorded_entry){
		table_row++;
		table_col = 0;
	}


	//step2 -> U(i,1) is known so fill up row(j)
	cout << "U(i,1) ------> \n\n";
	recorded_entry = 0;
	j = 1;
	for(i=0; i<=num_i; i++){
		if(i==0 || i==num_i)
			val = 0;
		else{
			if(constant_coeff){
				val = func_f(i*h) + k*func_g(i*h) + func_v_constant((i+1)*h)*func_v_constant((i+1)*h)*con*func_f((i+1)*h)/2.0 - func_v_constant(i*h)*func_v_constant(i*h)*con*func_f(i*h) + func_v_constant((i-1)*h)*func_v_constant((i-1)*h)*con*func_f((i-1)*h)/2.0;
			}
			else{
				val = func_f(i*h) + k*func_g(i*h) + func_v((i+1)*h)*func_v((i+1)*h)*con*func_f((i+1)*h)/2.0 - func_v(i*h)*func_v(i*h)*con*func_f(i*h) + func_v((i-1)*h)*func_v((i-1)*h)*con*func_f((i-1)*h)/2.0;
			}

			if(abs(val)<tol)
				val = 0;
		}

		//set row j to val at this i
		rowj.set(i,0, val);

		//possibly print table entry
		if( fmod((double)i,(double)delx)<tol && fmod((double)j,(double)delt)<tol ){
			Table.set(table_row,table_col, val);
			table_col++;
			recorded_entry = 1;
		}
	}
	//update table row and col
	if(recorded_entry){
		table_row++;
		table_col = 0;
	}


	//step3: iterate the 5 point scheme
	cout << "U(i,j+1) ------> \n\n";
	long vj;
	double p1,p2,p3,p4;
	for(j=1; j<num_j; j++){
		vj = j+1;
		recorded_entry = 0;
		for(i=0; i<=num_i; i++){

			if(i==0 || i==num_i){
				val = 0;
			}
			else{
				p1 = rowj.get(i,0);
				p2 = rowj.get(i+1,0);
				p3 = rowj.get(i-1,0);
				p4 = rowj1.get(i,0);

				if(constant_coeff){
					val = 2*( 1 - func_v_constant(i*h)*func_v_constant(i*h)*con )*p1 + func_v_constant((i+1)*h)*func_v_constant((i+1)*h)*con*p2 + func_v_constant((i-1)*h)*func_v_constant((i-1)*h)*con*p3 - p4;
				}
				else{
					val = 2*( 1 - func_v(i*h)*func_v(i*h)*con )*p1 + func_v((i+1)*h)*func_v((i+1)*h)*con*p2 + func_v((i-1)*h)*func_v((i-1)*h)*con*p3 - p4;
				}

				if(abs(val)<tol)
					val = 0;
			}

			//save current value U(i,j+1)
			row_temp.set(i,0, val);

			//possibly print table entry
			if( fmod((double)i,(double)delx)<tol && fmod((double)vj,(double)delt)<tol ){
				Table.set(table_row,table_col, val);
				table_col++;
				recorded_entry = 1;
			}
		}

		//update table row and col
		if(recorded_entry){
			table_row++;
			table_col = 0;
		}

		//at the end of i cycle for this j,
		//move up the two rows -> row(j-1) = row(j); row(j) = val = row(j+1)
		for(i=0; i<=num_i; i++){
			temp = rowj.get(i,0);
			val = row_temp.get(i,0);
			rowj1.set(i,0, temp);
			rowj.set(i,0, val);
		}

		//move on to next j..
	}

	//return the table ------------------>
	return Table;
}



/* returns a table of approximated values via finite differences (non-std algorithm) */
STLMatrix<double> finiteDiff::getNONSTDTable(double table_x_spacing_val, double table_t_spacing_val){

	table_x_spacing = table_x_spacing_val;
	table_t_spacing = table_t_spacing_val;

	//setup finite difference variables and storage containers
	long num_i = (long)ceil(L/h);
	long num_j = (long)ceil(T/k);

	cout << "setting up structures.. " << endl;

	//setup storage containers for the 3 rows we need
	STLMatrix<double> rowj(num_i+1,1);		//row j;
	STLMatrix<double> rowj1(num_i+1,1);		//row j-1;
	STLMatrix<double> row_temp(num_i+1,1);	//temp row

	//setup table variables
	long delx = (long)ceil(table_x_spacing/h);
	long delt = (long)ceil(table_t_spacing/k);
	long table_rows = (long)ceil(L/table_x_spacing);
	long table_cols = (long)ceil(T/table_t_spacing);
	long table_row = 0, table_col = 0;
	int recorded_entry = 0;
	STLMatrix<double> Table(table_rows+1,table_cols+1);
	

	cout << "starting finite differences.. " << endl;

	//perform finite differences -------------------->
	
	//step1 -> U(i,0) = f(ih) so fill up row(j-1)
	cout << "U(i,0) ------> \n\n";
	recorded_entry = 0;
	j = 0;
	for(i=0; i<=num_i; i++){
		if(i==0 || i==num_i)
			val = 0;
		else{
			val = func_f(i*h);
			if(abs(val)<tol)
				val = 0;
		}
		
		//set row j-1 to val at this i
		rowj1.set(i,0, val);

		//possibly print table entry
		if( fmod((double)i,(double)delx)<tol && fmod((double)j,(double)delt)<tol ){
			Table.set(table_row,table_col, val);
			table_col++;
			recorded_entry = 1;
		}
	}
	//update table row and col
	if(recorded_entry){
		table_row++;
		table_col = 0;
	}


	//step2 -> U(i,1) is known so fill up row(j)
	cout << "U(i,1) ------> \n\n";
	recorded_entry = 0;
	j = 1;
	for(i=0; i<=num_i; i++){
		if(i==0 || i==num_i)
			val = 0;
		else{
			if(constant_coeff){
				val = func_f(i*h) + k*func_g(i*h) + func_u_constant((i+1)*h)*func_u_constant((i+1)*h)*func_f((i+1)*h)/2.0 - func_u_constant(i*h)*func_u_constant(i*h)*func_f(i*h) + func_u_constant((i-1)*h)*func_u_constant((i-1)*h)*func_f((i-1)*h)/2.0;
			}
			else{
				val = func_f(i*h) + k*func_g(i*h) + func_u((i+1)*h)*func_u((i+1)*h)*func_f((i+1)*h)/2.0 - func_u(i*h)*func_u(i*h)*func_f(i*h) + func_u((i-1)*h)*func_u((i-1)*h)*func_f((i-1)*h)/2.0;
			}

			if(abs(val)<tol)
				val = 0;
		}

		//set row j to val at this i
		rowj.set(i,0, val);

		//possibly print table entry
		if( fmod((double)i,(double)delx)<tol && fmod((double)j,(double)delt)<tol ){
			Table.set(table_row,table_col, val);
			table_col++;
			recorded_entry = 1;
		}
	}
	//update table row and col
	if(recorded_entry){
		table_row++;
		table_col = 0;
	}


	//step3: iterate the 5 point scheme
	cout << "U(i,j+1) ------> \n\n";
	long vj;
	double p1,p2,p3,p4;
	for(j=1; j<num_j; j++){
		vj = j+1;
		recorded_entry = 0;
		for(i=0; i<=num_i; i++){

			if(i==0 || i==num_i){
				val = 0;
			}
			else{
				p1 = rowj.get(i,0);
				p2 = rowj.get(i+1,0);
				p3 = rowj.get(i-1,0);
				p4 = rowj1.get(i,0);

				if(constant_coeff){
					val = 2*( 1 - func_u_constant(i*h)*func_u_constant(i*h) )*p1 + func_u_constant((i+1)*h)*func_u_constant((i+1)*h)*p2 + func_u_constant((i-1)*h)*func_u_constant((i-1)*h)*p3 - p4;
				}
				else{
					val = 2*( 1 - func_u(i*h)*func_u(i*h) )*p1 + func_u((i+1)*h)*func_u((i+1)*h)*p2 + func_u((i-1)*h)*func_u((i-1)*h)*p3 - p4;
				}

				if(abs(val)<tol)
					val = 0;
			}

			//save current value U(i,j+1)
			row_temp.set(i,0, val);

			//possibly print table entry
			if( fmod((double)i,(double)delx)<tol && fmod((double)vj,(double)delt)<tol ){
				Table.set(table_row,table_col, val);
				table_col++;
				recorded_entry = 1;
			}
		}

		//update table row and col
		if(recorded_entry){
			table_row++;
			table_col = 0;
		}

		//at the end of i cycle for this j,
		//move up the two rows -> row(j-1) = row(j); row(j) = val = row(j+1)
		for(i=0; i<=num_i; i++){
			temp = rowj.get(i,0);
			val = row_temp.get(i,0);
			rowj1.set(i,0, temp);
			rowj.set(i,0, val);
		}

		//move on to next j..
	}

	//return the table ------------------>
	return Table;
}







/* f(x) <-> initial condition */
double finiteDiff::func_f(double x){
	double res;
	double pi = M_PI;
	double m = M;
	double l = 1.0;
	res = sin(m*pi*x);
	return res;
}


/* g(x) <-> initial condition */
double finiteDiff::func_g(double x){
	double res;
	res = 0;
	return res;
}


/* v(x) for non-exact (non-constant v) case */
double finiteDiff::func_v(double x){
	double res;
	res = 1.0+x*x;
	return res;
}


/* v(x) for exact (constant v) case */
double finiteDiff::func_v_constant(double x){
	double res;
	res = V;
	return res;
}


/* u(x) for non-exact (non-constant v) case */
double finiteDiff::func_u(double x){
	double res;
	double pi = M_PI;
	double m = M;
	double w = func_v(x)*m*pi;
	double d = w/func_v(x);

	res = sin((w*k)/2.0)/sin((d*h)/2.0);
	return res;
}


/* u(x) for exact (constant v) case */
double finiteDiff::func_u_constant(double x){
	double res;
	double pi = M_PI;
	double m = M;
	double w = func_v_constant(x)*m*pi;
	double d = w/func_v_constant(x);

	res = sin((w*k)/2.0)/sin((d*h)/2.0);
	return res;
}

