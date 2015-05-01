/*
 * Data.h
 *
 *  Created on: Apr 29, 2015
 *      Author: mwahab
 */

#ifndef DATA_H_
#define DATA_H_

#include <vector>

using namespace std;

class Point {

	double x;
	double y;

	public:

		Point(double x, double y){
			this->x = x;
			this->y = y;
		}

		double get_x(){
			return this->x;
		}

		double get_y(){
			return this->y;
		}
};

class Data {

	int index;
	int point_count;
	vector<Point> points;
	vector<double> xS;
	vector<double> yS;

public:
	Data();
	void add_point(double x, double y);
	Point get_point();
	Point get_next();
	vector<double> get_Xs();
	vector<double> get_Ys();
	int size();
	bool empty();
	virtual ~Data();
};

#endif /* DATA_H_ */
