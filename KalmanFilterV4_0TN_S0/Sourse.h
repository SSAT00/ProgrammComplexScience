#pragma once
#include"Kepler.h"

extern "C" __declspec(dllexport) int reform_start_data_and_run_epoch(double*, double*, double*, double*, double*);
extern "C" __declspec(dllexport) void reform_start_deta_and_make_start_observations(double*, double*, double*, double*);
