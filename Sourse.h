#pragma once
#include<math.h>
#include<iostream>

using namespace std;
extern "C" __declspec(dllexport) void Kepler2Pos(double*, double*);
extern "C" __declspec(dllexport) void Pos2Kepler(double*, double*);
