#ifndef ATTACH_H
#include "Attach.h"
#endif

/* 1D Problem */

// Solve y'(x) = 1 + (y/x) + (y/x)^{2}, y(x_{0}) = \alpha, x_{0} >= 1

// input is assumed to be of the form (t, x)

double examples::f1(double *x, int n)
{
	double t=(x[2]/x[1]); 

	return (1.0 + t + template_funcs::DSQR(t) ); 
}

/* 2D Problem */

// Solve y''(x) = \epsilon y'(x) ( 1 - ( y'(x) )^{2} ) - y, y(x_{0}) = \alpha, y'(x_{0}) = \beta, \epsilon > 0
// Break system into a pair of coupled first order equations
// Let u_{1} = y, u_{2} = y' => u_{1}'(t) = y'(t), u_{2}'(t) = y''(t)
// System to be solved is u_{1}'(t) = u_{2}(t), u_{2}'(t) = \epsilon u_{2}(t) ( 1 - ( u_{2}(t) )^{2} ) - u_{1}(t)

// input is assumed to be of the form (t, x, y)

double examples::g11(double *x, int n)
{
	return (x[3]); 
}

double examples::g12(double *x, int n)
{
	double epsilon = 0.1; 

	return ( epsilon * x[3] * (1.0 - template_funcs::DSQR(x[3]) ) - x[2] ); 
}

/* 3D Problem */
// Solve the system of three coupled first order equations
// x'(t) = y(t) z(t), y'(t) = -2 x(t) z(t), z'(t) = x(t) y(t), 
// x(t_{0}) = \alpha, y(t_{0}) = \beta, z(t_{0}) = \gamma, t_{0} >= 0

// input is assumed to be of the form (t, x, y, z)

double examples::h11(double *x, int n)
{
	return (x[3]*x[4]); 
}

double examples::h12(double *x, int n)
{
	return (-2.0*x[2]*x[4]); 
}

double examples::h13(double *x, int n)
{
	return (x[2]*x[3]); 
}

void examples::solve_p1()
{
	int nn; 

	double a, b, ic; 

	/**********************************************************/
	/* Problem 1 */

	// Solve y'(x) = 1 + (y/x) + (y/x)^{2}, y(x_{0}) = \alpha, x_{0} >= 1
	
	nn = 21; a=1.0; b=4.0; ic=0.0; 

	IVP_Solver test(nn,a,b,ic,f1);

	/*test.RK4_Solve(); 
	test.output_solution("Pr1_Sol_RK4.txt");*/

	test.ABM_PC_Solve();
	test.output_solution("Pr1_Sol_PC.txt");
}

void examples::solve_p2()
{
	/**********************************************************/
	/* Problem 2 */

	// Solve y''(x) = \epsilon y'(x) ( 1 - ( y'(x) )^{2} ) - y, y(x_{0}) = \alpha, y'(x_{0}) = \beta, \epsilon > 0
	// Break system into a pair of coupled first order equations
	// Let u_{1} = y, u_{2} = y' => u_{1}'(t) = y'(t), u_{2}'(t) = y''(t)
	// System to be solved is u_{1}'(t) = u_{2}(t), u_{2}'(t) = \epsilon u_{2}(t) ( 1 - ( u_{2}(t) )^{2} ) - u_{1}(t)
	// The nice thing about solving second order odes numerically is that you get y'(t) for free when you compute y(t)

	int nn; 

	double a, b, ic1, ic2; 

	nn=301; a=0.0; b = 50.0; ic1=3.0; ic2=2.0; 

	IVP_Solver test1(nn,a,b,ic1,ic2,g11,g12);

	/*test1.RK4_Solve(); 
	test1.output_solution("Pr2_Sol_RK4.txt");*/

	test1.ABM_PC_Solve(); 
	test1.output_solution("Pr2_Sol_PC.txt");
}

void examples::solve_p3()
{
	int nn; 

	double a, b, ic1, ic2, ic3; 

	/**********************************************************/
	/* Problem 3 */
	// Solve the system of three coupled first order equations
	// x'(t) = y(t) z(t), y'(t) = -2 x(t) z(t), z'(t) = x(t) y(t), 
	// x(t_{0}) = \alpha, y(t_{0}) = \beta, z(t_{0}) = \gamma, t_{0} >= 0

	nn = 301; a = 0.0; b = 50.0; ic1 = 2.0; ic2 = 1.0; ic3 = 0.5; 

	IVP_Solver test2(nn, a, b, ic1, ic2, ic3, h11, h12, h13);

	/*test2.RK4_Solve(); 
	test2.output_solution("Pr3_Sol_RK4.txt");*/

	test2.ABM_PC_Solve(); 
	test2.output_solution("Pr3_Sol_PC.txt");
}