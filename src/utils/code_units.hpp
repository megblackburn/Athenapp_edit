#ifndef CODE_UNITS_H
#define CODE_UNITS_H

#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <iomanip>

#include "../athena.hpp"
#include "../athena_arrays.hpp"

static const Real CONST_pc = 3.086e18;
static const Real CONST_yr = 3.154e7;
static const Real CONST_amu = 1.66053886e-24;
static const Real CONST_kB = 1.3806505e-16;
static const Real unit_length = CONST_pc*1e3;
static const Real unit_time = CONST_yr*1e6;
static const Real unit_density = CONST_amu;
static const Real unit_velocity = unit_length/unit_time;
static const Real unit_q = (unit_density * pow(unit_velocity, 3.0))/unit_length;
static const Real KELVIN = unit_velocity*unit_velocity*CONST_amu/CONST_kB;

static Real Zsol = 1.0;
static Real Xsol = 1.0;

static Real X = Xsol * 0.7381;
static Real Z = Zsol * 0.0134;
static Real Y = 1 - X - Z;

static Real mu = 1.0/(2.0*X+3.0*(1.0-X-Z)/4.0+Z/2.0);
static Real mue = 2.0/(1.0+X);
static Real muH = 1.0/X;

static Real T_floor = 1e4;
static Real T_ceil = 1e8;
static Real T_hot = 1e6;
static Real T_hot_req = 1e6;
static Real T_cold = 1e4;
static Real T_cut_mul = 0.6;
static Real T_cut = T_cut_mul*T_hot_req;


#endif
