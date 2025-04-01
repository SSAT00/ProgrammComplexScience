#pragma once
#include"Eigen/Dense"
#include<iostream>
#include<iomanip>
#include<map>
#include<random>
#include<thread>
#include<fstream>

using namespace std;
using namespace Eigen;

void KeplerX(Matrix<double, 6, 1>, double, Matrix<double, 6, 1>*);
void KeplerXDX(Matrix<double, 6, 1>, double, Matrix<double, 6, 1>*, Matrix<double, 6, 6>*);
