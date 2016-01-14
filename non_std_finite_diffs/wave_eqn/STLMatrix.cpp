/* 
 * Simple arbitrary data type matrix class using the Standard Template Library
 *
 * Sergey Voronin */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>


using namespace std;

template <typename T> class STLMatrix {

	private:

		vector< vector<T> > MatVect; //a vector of vectors of data type T representing a matrix
		long num_rows;
		long num_cols;

	public:
	
		/* constructors */
		STLMatrix();
		STLMatrix(long nrows, long ncols);
	
		/* utility functions */
		void resizeMatrix(long nrows, long ncols);
		void set(long row, long col, T element);
		long getNumRows();
		long getNumCols();
		T get(long row, long col);
		T max();
		double mean();
		STLMatrix<T> abs();
		STLMatrix<T> transpose();
	
		/* overloaded operators */
                friend STLMatrix<T> operator+<>(const STLMatrix<T> &A, const STLMatrix<T> &B);
		friend STLMatrix<T> operator-<>(const STLMatrix<T> &A, const STLMatrix<T> &B);
		//friend STLMatrix<T> operator*<>(const T constant, const STLMatrix<T> &A);
		friend ostream& operator<<<>(ostream& s, const STLMatrix<T> &A);

};


/* empty constructor */
template <typename T> STLMatrix<T>::STLMatrix(){
	resizeMatrix(1, 1);
	num_rows = 1;
	num_cols = 1;
}


/* constructor with parameters */
template <typename T> STLMatrix<T>::STLMatrix(long nrows, long ncols){
	resizeMatrix(nrows, ncols);
	num_rows = nrows;
	num_cols = ncols;
}


/* resizes the MatVect matrix */
template <typename T> void STLMatrix<T>::resizeMatrix(long nrows, long ncols){
	MatVect.resize(nrows, vector<T>(ncols));
	num_rows = nrows;
	num_cols = ncols;
}


/* returns the number of rows of the matrix */
template <typename T> long STLMatrix<T>::getNumRows(){
	return num_rows;
}


/* returns the number of columns of the matrix */
template <typename T> long STLMatrix<T>::getNumCols(){
	return num_cols;
}


/* sets the value of a specific element */
template <typename T> void STLMatrix<T>::set(long row, long col, T element){
	if((row>=0 && row<num_rows) && (col>=0 && col<num_cols))
		MatVect[row][col] = element;
	else{
		cout << "Out of range!\n" << endl;
		cout << row << "\t" << col << endl;
	}
}


/* returns the element of the matrix at the specified position */
template <typename T> T STLMatrix<T>::get(long row, long col){
	T val;
	if((row>=0 && row<num_rows) && (col>=0 && col<num_cols))
		val =  MatVect[row][col];
	else{
		cout << "Out of range!\n" << endl;
		cout << row << "\t" << col << endl;
		val = 0;
	}
	return val;
}


/* returns the maximum value of the matrix */
template <typename T> T STLMatrix<T>::max(){
	T val, max;
	max = MatVect[0][0];
	for(long row=0; row<num_rows; row++){
        for(long col=0; col<num_cols; col++){
			val = MatVect[row][col];
			if(val >= max)
				max = val;
		}
	}
	return max;
}


/* returns the mean value of the matrix */
template <typename T> double STLMatrix<T>::mean(){
	T val, sum;
	double mean;
	sum = 0;
	for(long row=0; row<num_rows; row++){
        for(long col=0; col<num_cols; col++){
			val = MatVect[row][col];
			sum += val;
		}
	}
	mean = sum/(num_rows*num_cols);
	return mean;
}


/* returns a new matrix B where each b = abs(a) of each a of matrix A */
template <typename T> STLMatrix<T> STLMatrix<T>::abs(){
	STLMatrix<T> matA(num_cols, num_rows);
	T val;
	for(long row=0; row<num_rows; row++){
        for(long col=0; col<num_cols; col++){
			val = MatVect[row][col];
			if(val < (T)0)
				matA.MatVect[row][col] = (T)(-1) * MatVect[row][col];
			else
				matA.MatVect[row][col] = MatVect[row][col];
		}
	}
	return matA;
}


/* returns the transpose of the matrix */
template <typename T> STLMatrix<T> STLMatrix<T>::transpose(){
	STLMatrix<T> matTA(num_cols, num_rows); //holds the transpose matrix TA
	for(long row=0; row<num_rows; row++){
        for(long col=0; col<num_cols; col++){
			matTA.MatVect[col][row] = MatVect[row][col];
		}
	}
	return matTA;
}



/* addition operator - performs the addition of two matrices A and B and returns a new matrix C = A + B */
template <typename T> STLMatrix<T> operator+(const STLMatrix<T> &A, const STLMatrix<T> &B){
	STLMatrix<T> matC(A.num_rows, A.num_cols);
	int row, col;
	
	for(row=0; row<A.num_rows; row++){
		for(col=0; col<A.num_cols; col++){
			matC.MatVect[row][col] = A.MatVect[row][col] + B.MatVect[row][col];	
		}
	}
	
	return matC;
}


/* subtraction operator - performs the subtraction of matrix B from A and returns a new matrix C = A - B */
template <typename T> STLMatrix<T> operator-(const STLMatrix<T> &A, const STLMatrix<T> &B){
	STLMatrix<T> matC(A.num_rows, A.num_cols);
	int row, col;
	
	for(row=0; row<A.num_rows; row++){
		for(col=0; col<A.num_cols; col++){
			matC.MatVect[row][col] = A.MatVect[row][col] - B.MatVect[row][col];	
		}
	}
	
	return matC;
}


/* multiplication by a constant operator - returns new matrix B = constant*A */
/*
template <typename T> STLMatrix<T> operator*(T constant, const STLMatrix<T> &A){
	STLMatrix<T> matB(A.num_rows, A.num_cols);
	int row, col;
	
	for(row=0; row<A.num_rows; row++){
		for(col=0; col<A.num_cols; col++){
			matB.MatVect[row][col] = A.MatVect[row][col] * constant;	
		}
	}
	
	return matB;
}*/


/* prints out the current matrix to the ouput stream*/
template <typename T> ostream& operator<<(ostream& s, const STLMatrix<T> &A){
	long row, col, endcol;
	char *separator = " ";
	
	if(A.num_cols>1){
		for(row=0; row<A.num_rows; row++){
			for(col=0; col<A.num_cols; col++){
				s << A.MatVect[row][col] << separator;
			}
			cout << endl;
		}
	}
	
	else if(A.num_cols==1){
		for(row=0; row<A.num_rows; row++)
			s << A.MatVect[row][0] << endl;
		s << endl;
	}
	s << endl;
}
