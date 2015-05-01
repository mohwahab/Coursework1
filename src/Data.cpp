/*
 * Data.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: mwahab
 */

#include "Data.h"

Data::Data() {
	// TODO Auto-generated constructor stub
	index = 0;
	point_count = 0;
}

Data::~Data() {
	// TODO Auto-generated destructor stub
}

bool Data::empty(){
	return this->points.empty();
}

void Data::add_point(double x, double y){
	Point point(x, y);
	this->xS.push_back(x);
	this->yS.push_back(y);
	this->points.push_back(point);
}

Point Data::get_point(){
	Point point = this->points.back();
	this->points.pop_back();
	return point;
}

Point Data::get_next(){

	Point point(xS[index],yS[index]);
	index++;
	if(index >= point_count){
		index = 0;
	}
	return point;
}

vector<double> Data::get_Xs(){
	return this->xS;
}

vector<double> Data::get_Ys(){
	return this->yS;
}

int Data::size(){
	return point_count;
}
