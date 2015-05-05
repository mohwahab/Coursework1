//============================================================================
// Name        : Regression.cpp
// Author      : M.Wahab
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdlib.h> 	// std::atof
#include <cmath>        // std::abs
#include <algorithm>	// std::count
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Data.h"
#include "util/Matrix.h"
#include "gnuplot_i.hpp"

using namespace std;

Matrix<double>* load_data(char* filename, vector<double>& b){


	int delimiter_index = 0;
	int row = 0, col = 0;
	int size = -2;
	char delimiter = ',';
	double row_elem;

	ifstream file (filename);
	string line, sub_str, a_row, b_row;
	Matrix<double>* a_ptr;
	cout << "***** Going to load Matrix data from file ["<<filename<<"] *****"<< endl; // prints Hello World!!!

		while (getline(file, line))
		{

			 delimiter_index = line.find(delimiter);
			 a_row = line.substr(0, delimiter_index);
			 b_row = line.substr(delimiter_index+1);


			 if(size == -2){

				 size = count(a_row.begin(), a_row.end(), ' ')+1;
				 a_ptr = new Matrix<double>(size, size);

			 }

			 b.push_back(atof(b_row.c_str()));

			 istringstream a_row_ss(a_row);

			 col = 0;

			 while(a_row_ss >> row_elem){

				 (*a_ptr)(row, col) =  row_elem;
				 col++;
			 }

			 row++;

		}
		cout<<"ROW COUNT = "<<(a_ptr->get_rows())<<" COL COUNT = "<<(a_ptr->get_cols())<<endl;
		file.close();
		return a_ptr;
}

Data load_data(char* filename){

	double x, y;
	ifstream file (filename);
	string line;
	//string input;
	Data data;

	cout << "***** Going to load data from file ["<<filename<<"] *****"<< endl; // prints Hello World!!!

	while (getline(file, line))
	{
	     istringstream ss(line);
	     /*
	     getline(ss, input, ' ');
	     x = atof(input.c_str());
	     getline(ss, input, ' ');
	     y = atof(input.c_str());
	     */

	     if (!(ss >> x >> y)) { // error
	    	 cerr << "Invalid data format: [" << line << "]" <<endl;
	    	 break;
	     }

	     data.add_point(x, y);
	     cout << "[X= "<< x <<", Y= " << y <<"]"<<endl;
	}

	file.close();

	return data;
}

void linear_regression(Data data, int presist_for){

	vector<double> points_x;
	vector<double> points_y;

	double 	sum_x = 0,
			sum_y = 0,
			sum_xy = 0,
			sum_x2 = 0;
	double 	st = 0,
			sr = 0;
	double x, y;

	//int n = 0; //Number of points

	cout << "***** Going to process data (Linear Regression) *****" << endl; // prints Hello World!!!

	while(!data.empty()){

		Point point = data.get_point();

		x = point.get_x();
		y = point.get_y();

		points_x.push_back(x);
		points_y.push_back(y);

		cout << "[X= " << point.get_x() << ", Y= " << point.get_y() <<"]"<< endl;

		sum_x += x;
		sum_y += y;
		sum_xy += x*y,
		sum_x2 += x*x;

		//n++;
	}

	int n = points_x.size(); //Total number of points

	double xm = sum_x/n; //Mean values of x
	double ym = sum_y/n; //Mean values of y

	double a1 = (n*sum_xy - sum_x*sum_y)/(n*sum_x2 - sum_x*sum_x);
	double a0 = ym - a1*xm;

	double 	yi_ym = 0,
			yi_a0_a1xi = 0;

	for (int i= 0; i < n; i++) {

		yi_ym = points_y[i] - ym;
		yi_a0_a1xi = points_y[i] - a0 - a1*points_x[i] ;

		st += yi_ym*yi_ym;
		sr += yi_a0_a1xi*yi_a0_a1xi;
	}

	cout << "[a0=" << a0 << "]" << "[a1=" << a1 << "]" << endl;
	cout << "[St=" << st << "]" << "[Sr=" << sr << "]" << endl;

	//Gnuplot gnu_plot = Gnuplot("Linear Regression", "points", "X-axis", "Y-axis", points_x, points_y);

	std::ostringstream a0_strs;
	a0_strs << a0;
	std::string a0_str = a0_strs.str();

	std::ostringstream a1_strs;
	a1_strs << a1;
	std::string a1_str = a1_strs.str();

	string title = "Linear Regression";

	Gnuplot gnu_plot = Gnuplot("lines");

	gnu_plot.set_style("points");
	gnu_plot.plot_xy(points_x, points_y, title);

	gnu_plot.set_style("lines");
	gnu_plot.set_xlabel("Y-axis");
	gnu_plot.set_ylabel("X-axis");
	gnu_plot.plot_slope(a1, a0, "[ y = "+a0_str+" + "+a1_str+" x ]");

	/*
	Gnuplot g1 = Gnuplot("lines");

	cout << endl << endl << "*** user-defined lists of points" << endl;

	 vector<double> xx;
	 vector<double> yy;

	for (int i = 0; i < 50; i++)
	  {
		xx.push_back((double)i);
		yy.push_back((double)i * (double)i);
	  }

	g1.reset_plot();
	g1.set_style("points");
	g1.plot_xy(xx,yy,"user-defined points");
	*/

	sleep(presist_for);
}


void newtons_interpolation(Data data, double xi){

	int n = data.size();
	double fdd[n][n];


	vector<double> x = data.get_Xs();
	vector<double> y = data.get_Ys();

	//Fill 1st column of the table --> f(x0), f(x1), ..... f(xn-1)
	for (int i = 0; i < n; i++) {
		fdd[i][0] = y[i];
		cout<<"fdd["<<i<<"][0] = "<<fdd[i][0]<<endl;
	}

	//Calculate the divided differences and store the  values in table --> f[xi,xj] = f(xi)-f(xj)/xi-xj  &  f[xi,xj,xk] = f[xi,xj] - f[xj,xk] / xi-xk
	for (int j = 1; j < n; j++) {

		for (int i = 0; i < n-j; i++) {
			fdd[i][j] = (fdd[i+1][j-1] - fdd[i][j-1]) / (x[i+j] - x[i]); //f[x1,x0] = f(x1) - f(x0)/x1-x0  --- 2D-Array Notation ---> f[0][1] = f[1][0] - f[0][0] / x[1]-x[0]
			cout<<"fdd["<<i<<"]["<<j<<"] = "<<fdd[i][j]<<endl;
		}

	}
	//After this loop we should have the table filled so that we can calculate all the coefficients

	double xterm = 1; 	//Polynomial terms (x-x0)(x-x1) ..... (x-xn-1)
	double yout = 0; 	//To hold the output value
	double yint[n];   	//Array to hold the terms of polynomial equation
	double ea[n];		//Error array, to hold error in each term

	yint[0] = fdd[0][0]; //y(x0) = f(x0)

	std::ostringstream yint_stream;
	yint_stream << "" << yint[0];

	string equation = yint_stream.str();
	//string xterm_str("");

	std::ostringstream xterm_stream, equation_stream;

	for (int order = 1; order < n; order++) {
		xterm = xterm * (xi - x[order-1]);
		xterm_stream << "*(x-" << x[order-1] << ")";
		equation_stream << "+" << fdd[0][order] << xterm_stream.str();
		yout = yint[order-1] + fdd[0][order] * xterm;
		ea[order-1] = yout - yint[order-1];
		yint[order] = yout;
		//cout<<"yint["<<order<<"] = "<<yint[order]<<endl;
	}

	equation += equation_stream.str();
	cout << "Y[x] = " << equation << endl;
	cout << "Y["<<xi<<"] = " << yout << endl;

	string title = "Newton's Interpolation";
	Gnuplot gnu_plot = Gnuplot("lines");
	gnu_plot.set_style("points");
	gnu_plot.plot_xy(x, y, title);

	gnu_plot.set_style("lines");
	gnu_plot.set_xlabel("Y-axis");
	gnu_plot.set_ylabel("X-axis");
	gnu_plot.plot_equation(equation, equation/*title*/);

	sleep(30);
}


/************************** Gauss Elemination (Page 264)*****************************/

void pivot (Matrix<double>& a, vector<double>& b, vector<double>& s, int n, int k){

	int p = k;
	double dummy;
	double big = abs(a(k,k)/s[k]);

	for (int ii = k+1; ii <= n; ii++) {

		dummy = abs(a(ii,k)/s[ii]);

		if(dummy > big){
			big = dummy;
			p = ii;
		}
	}

	if(p != k){

		for( int jj = k; jj <= n; jj++){
			dummy = a(p,jj);
			a(p,jj) = a(k,jj);
			a(k,jj) = dummy;
		}

		dummy = b[p];
		b[p] = b[k];
		b[k]= dummy;
		dummy = s[p];
		s[p] = s[k];
		s[k] = dummy;
	}
}

int eliminate(Matrix<double>& a, vector<double>& s, int n, vector<double>& b, double tol){

	int er = 0;

	double factor = 1;



	for (int k = 0; k <= n; k++) {

		pivot(a, b, s, n, k);
		cout<<"********* PIVOT["<<k<<"]**********"<<endl;
		for (int i=0; i<a.get_rows(); i++) {
				for (int j=0; j<a.get_cols(); j++) {
				  cout << a(i,j) << " ";
				}
				cout<<","<<b[i]<<endl;;
			  }

		cout<<"*******************"<<endl;

		if(abs(a(k,k)/s[k]) < tol){
			er = -1;
			return er;
		}

		for (int i = k+1; i <= n; i++) {

				factor = a(i,k)/a(k,k);

				for (int j = k; j <= n; j++) {
					cout<<"HERE3-7"<<endl;
					a(i,j) = a(i,j) - factor*a(k,j);
					cout<<"HERE3-8 j= "<<j<<endl;
				}

				b[i] = b[i] - factor * b[k];
				cout<<">>>>>>>>>>>>>>>> HERE3-9 i= "<<i<<" k= "<<k<<endl;
		}
	}


	if(abs(a(n,n)/s[n]) < tol){
		er = -1;
	}

	return er;
}

/*vector<double>*/ void substitute (Matrix<double> a, int n, vector<double> b,vector<double>& x){

	cout<<"*******************"<<endl;
	for (int i=0; i<a.get_rows(); i++) {
		    for (int j=0; j<a.get_cols(); j++) {
		      cout << a(i,j) << " ";
		    }
		    cout<<","<<b[i]<<endl;;
		  }

	cout<<"*******************"<<endl;
	int sum = 0;

	//vector<double> x;
	cout<<"HERE4-1 >>> b["<<n<<"] = "<<b[n]<<endl;
	cout<<"HERE4-1 >>> a["<<n<<"] = "<<a(n,n)<<endl;
	cout<<"HERE4-1 >>> b/a = "<<a(n,n)<<endl;
	//x.push_back(b[n]/a(n,n));
	x[n] = b[n]/a(n,n);
	cout<<"HERE4-1 >>> x["<<n<<"] = "<<x[n]<<endl;
	for(int i = n - 1; i >= 1; --i){

		sum = 0;
		for(int j = i+1; j <= n; ++j){

			sum += a(i,j) * x[j];
			cout<<"HERE4-4 >>> x["<<j<<"] = "<<x[j]<<endl;
			cout<<"HERE4-4 j= "<<j<<endl;
		}
		cout<<"HERE4-5"<<endl;
		//x[n] = (b[n] - sum) / a(n,n);
		x[i] = (b[i] - sum) / a(i,i);

		cout<<">>>>>>>>>>>>>>>>>>>> HERE4-6 i= "<<i<<endl;
	}
	//x[0] = (b[3] - sum) / a(0,0);
	cout<<"SUBISTIUTE DONE"<<endl;
	//return x;
}



/*vector<double>*/ void gauss_elemination(Matrix<double> a, vector<double> b, int n, vector<double>& x, double tol = 0){

	vector<double> s(n+1,0); //Scale vector
	//vector<double> x;
	cout<<"HERE1"<<endl;
	for (int i = 0; i <= n; i++) {

		s[i] = abs(a(i,1));

		for (int j = 1; j <= n; j++) {
			cout<<"----------- i = "<<i<<" j = "<<j<<" a = "<<a(i,j)<<" s = "<<s[i]<<endl;
			if(abs(a(i,j)) > s[i]) {

				s[i] = abs(a(i,j));

			}
		}
	}

	for(int k=0; k< s.size(); k++){
		cout<<">>>>>> s= "<<s[k]<<endl;
	}


	int er = eliminate(a, s, n, b, tol); //if er == -1 --> Singular matrix has been detected

	if( er != -1){
		/*x = */substitute(a, n, b, x);

	}else{
		cerr<<"Singular matrix has been detected, will exit"<<endl;
	}
	for(int k=0; k<= n; k++){
			cout<<">>>>>> X = "<<x[k]<<endl;
		}
	//return x;
}

/**********************************************************************/

int main(int argc, char *argv[]) {

	//Data data = load_data(argv[1]);
	//linear_regression(data, atoi(argv[2]));
	//newtons_interpolation(data, atof(argv[2]));

	vector<double> B;
	Matrix<double>* A = load_data(argv[1], B);

	for (int i=0; i<A->get_rows(); i++) {
	    for (int j=0; j<A->get_cols(); j++) {
	      cout << (*A)(i,j) << " ";
	    }
	    cout <<","<<B[i]<< endl;
	  }

	vector<double> x(B.size(), 0);
	gauss_elemination((*A), B, B.size()-1, x);
	for(int k=0; k< B.size(); k++){
		cout<<"X"<<k<<" = "<<x[k]<<endl;
	}
	delete A;

	return 0;
}
