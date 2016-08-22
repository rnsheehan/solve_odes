#ifndef EXAMPLES_H
#define EXAMPLES_H

// Test functions for the IVP Solver

namespace examples{

	/* 1D Problem */
	/* input is assumed to be of the form (t, x) */

	double f1(double *x, int n); 

	/* 2D Problem */
	/* input is assumed to be of the form (t, x, y) */
	double g11(double *x, int n); 
	double g12(double *x, int n); 

	/* 3D Problem */
	/* input is assumed to be of the form (t, x, y, z) */
	double h11(double *x, int n); 
	double h12(double *x, int n); 
	double h13(double *x, int n); 

	void solve_p1(); 
	void solve_p2(); 
	void solve_p3(); 

}

#endif