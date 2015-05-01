//============================================================================
// Name        : Regression.cpp
// Author      : M.Wahab
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdlib.h> //For atof
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Data.h"
#include "gnuplot_i.hpp"

using namespace std;

Data load_data(char* filename){

	double x, y;
	ifstream file (filename);
	string line;
	string input;
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
	    	 cout << "Invalid data format: [" << line << "]" <<endl;
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

/*
void newtons_interpolation(Data data, double xi){

	int n = data.size();
	double fdd[n][n];


	vector<double> x = data.get_Xs();
	vector<double> y = data.get_Ys();

	for (int i = 0; i < n; i++) {
		fdd[i][0] = y[i];
	}

	for (int var = 0; var < max; ++var) {

	}
}
*/

int main(int argc, char *argv[]) {

	Data data = load_data(argv[1]);
	linear_regression(data, atoi(argv[2]));

	return 0;
}
