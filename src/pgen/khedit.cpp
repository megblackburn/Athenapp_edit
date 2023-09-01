//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kh.cpp
//! \brief Problem generator for KH instability.
//!
//! Sets up several different problems:
//!   - iprob=1: slip surface with random perturbations
//!   - iprob=2: tanh profile, with single-mode perturbation (Frank et al. 1996)
//!   - iprob=3: tanh profiles for v and d, SR test problem in Beckwith & Stone (2011)
//!   - iprob=4: tanh profiles for v and d, "Lecoanet" test
//!   - iprob=5: two resolved slip-surfaces with m=2 perturbation for the AMR test

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>      // log
#include <cstring>    // strcmp()
#include <fstream>
#include <iomanip>

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../defs.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"
#include "../utils/townsend_cooling_max.hpp"
#include "../utils/hst_func.hpp"
#include "../utils/code_units.hpp"

namespace {
Real vflow;
int iprob;
Real PassiveDyeEntropy(MeshBlock *pmb, int iout);
} // namespace

Real threshold;
int RefinementCondition(MeshBlock *pmb);

static Real tfloor, tnotcool, tcut_hst, r_drop;
static Real Lambda_fac, Lambda_fac_time;         // for boosting cooling
static Real total_cooling=0.0;

// Returns unique pointer
// This is a function in C++13 onwards
// The implementation here is copied from https://stackoverflow.com/a/17903225/1834164

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

std::unique_ptr<Cooling> cooler;
/* ------------------------------*/


// Variable read from the input file

static Real L1 = 0.0;
static Real L2 = 0.0;
static Real L3 = 0.0;

static int nx1 = 1;
static int nx2 = 1;
static int nx3 = 1;

static Real x3min = 0.0;
static Real x3max = 0.0;

static int cooling_flag = 1;
static int shift_flag   = 1;

static Real amb_rho = 1.0;

static Real front_thick = 2.5;
static Real v_shear = 100 * (1.023*1e-3); // 100 km/s in kpc/Myr
static Real v_shift0 = 0.0;                // velocity of shift in z

static Real front_velocity = 0.0;

static Real knx_KH  = 1.0;
static Real kny_KH  = 1.0;
static Real amp_KH = 0.01;    // Amplitude = amp_KH * v_shear

static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;

static Real front_posn_old = 0.0;
// static Real v_shift_t_old  = 0.0;

static bool cooling_flag_print_count = false;
static bool shift_flag_print_count   = false;

static bool DEBUG_FLAG_MIX = false;

static Real time_last  = -1.0;
static Real shift_start = 0.0;

static Real debug_hst = 0.0;

void read_input (ParameterInput *pin){
  L1 = pin->GetReal("mesh","x1max") - pin->GetReal("mesh","x1min");
  L2 = pin->GetReal("mesh","x2max") - pin->GetReal("mesh","x2min");
  L3 = pin->GetReal("mesh","x3max") - pin->GetReal("mesh","x3min");

  x3min = pin->GetReal("mesh","x3min");
  x3max = pin->GetReal("mesh","x3max");

  nx1 = pin->GetInteger("mesh","nx1");
  nx2 = pin->GetInteger("mesh","nx3");
  nx3 = pin->GetInteger("mesh","nx3");

  cooling_flag = pin->GetInteger("problem","cooling_flag");
  shift_flag   = pin->GetInteger("problem","shift_flag");

  shift_start  = pin->GetReal("problem","shift_start");


  amb_rho      = pin->GetReal("problem","amb_rho");

  front_thick  = pin->GetReal("problem","front_thickness");
  v_shear      = pin->GetReal("problem","v_shear");
  v_shift0      = pin->GetReal("problem","v_shift");

  knx_KH = pin->GetReal("problem","knx_KH");
  kny_KH = pin->GetReal("problem","kny_KH");
  amp_KH = pin->GetReal("problem","amp_KH");

  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_cold       = pin->GetReal("problem","T_cold");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");

  rho_cold = pin->GetReal("problem", "rho_cold");
  rho_hot = pin->GetReal("problem", "rho_hot");

  T_cut = T_cut_mul*T_hot;

  Xsol            = pin->GetReal("problem","Xsol");
  Zsol            = pin->GetReal("problem","Zsol");

  X = Xsol * 0.7381;
  Z = Zsol * 0.0134;
  Y = 1 - X - Z;

  mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
  mue = 2.0/(1.0+X);
  muH = 1.0/X;

  Lambda_fac     = pin->GetReal("problem","Lambda_fac");
  DEBUG_FLAG_MIX = pin->GetBoolean("problem","DEBUG_FLAG");

  if (MAGNETIC_FIELDS_ENABLED) {
    B_x = pin->GetReal("problem", "B_x");
    B_y = pin->GetReal("problem", "B_y");
    B_z = pin->GetReal("problem", "B_z");
  }

  printf("____ Input file read! ______\n");

  return;

}

/*---------------------------------- townsend cooling ------------------------*/
void townsend_cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  Real g = pmb->peos->GetGamma();

  // Real sim_time = pmb->pmy_mesh->time;

  int il = pmb->is - NGHOST; int iu = pmb->ie + NGHOST;
  int jl = pmb->js - NGHOST; int ju = pmb->je + NGHOST;
  int kl = pmb->ks - NGHOST; int ku = pmb->ke + NGHOST;

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {

        debug_hst += 1.0;

        Real temp = (prim(IPR,k,j,i) / prim(IDN,k,j,i)) * KELVIN * mu ;

        //TODO: Move NaN checks to Cooling header file
        if (std::isnan(temp)) {

          printf("[cooling] nan detected for T at i,j,k : (%d, %d, %d) \n", i,j,k);
          printf("temp = %lf \n", temp );
          printf("rho  = %lf \n", cons(IDN,k,j,i));
          printf("IM1  = %lf \n", cons(IM1,k,j,i));
          printf("IM2  = %lf \n", cons(IM2,k,j,i));
          printf("IM3  = %lf \n", cons(IM3,k,j,i));

          Real E_kin  = SQR(prim(IVX,k,j,i));
          E_kin      += SQR(prim(IVY,k,j,i));
          E_kin      += SQR(prim(IVZ,k,j,i));
          E_kin      *= 0.5*prim(IDN,k,j,i);

          Real E_mag = 0.0;

          if (MAGNETIC_FIELDS_ENABLED) {
            E_mag += prim(IB1,k,j,i)*prim(IB1,k,j,i);
            E_mag += prim(IB2,k,j,i)*prim(IB2,k,j,i);
            E_mag += prim(IB3,k,j,i)*prim(IB3,k,j,i);
            E_mag *= 0.5;
          }

          // Set current temperature to T_floor
          cons(IEN,k,j,i) = E_kin + E_mag + (T_floor/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);
          continue;

        } // End of NaN check on current temperature

        if (temp > T_floor) {

          //* For own townsend cooling

          // Real *lam_para = Lam_file_read(temp);
          //printf("%f %f %f\n",lam_para[0],lam_para[1],lam_para[2]);
          // Real temp_new = T_new(temp, lam_para[0], lam_para[1], lam_para[2],
          //                         prim(IDN,k,j,i), dt, 
          //                         mu, mue, muH, 
          //                         g, T_floor, T_ceil, T_cut);
          // Defined in ../utils/townsend_cooling.hpp
          //  NOTE: Above, dt is in code units to avoid overflow. unit_time is cancelled in 
          //  calculation of T_new as we calculate (dt/tcool)
          // cons(IEN,k,j,i) += ((temp_new-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);

          //* For Max's townsend cooling

          if (temp<T_cut){
            Real rho      = cons(IDN,k,j,i);
            Real temp_cgs = temp;
            Real rho_cgs  = rho * unit_density;
            Real dt_cgs   = dt  * unit_time;
            Real cLfac    = Lambda_fac;

            Real temp_new = fmax(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, cLfac), T_floor);

            if  (std::isnan(temp_new)) {

              printf("[cooling] nan detected for T_new at i,j,k : (%d, %d, %d) \n", i,j,k);
              printf("temp = %lf \n", temp_cgs );
              printf("rho  = %lf \n", rho);
              printf("IM1  = %lf \n", cons(IM1,k,j,i));
              printf("IM2  = %lf \n", cons(IM2,k,j,i));
              printf("IM3  = %lf \n", cons(IM3,k,j,i));

              // temp_new set to T_floor
              Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);
              cons(IEN,k,j,i) += ccool;

            } //* End of NaN check

            else {

              Real ccool = ((temp_new-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);

              //TODO: Check the volume check with homogenous box
              //TODO: Put a test with ccool = 1

              cons(IEN,k,j,i) += ccool;

              if (DEBUG_FLAG_MIX){
                total_cooling -= 1.0;
              }
              else{
                total_cooling -= ccool;
              }


            } // End of else for NaN check

          } // End of T<T_cut check
        } // End of T>T_floor check

        else{  // If T<=T_floor

          // printf("+");

          Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);
          cons(IEN,k,j,i) += ccool;
          // total_cooling -= ccool;

        } // End of else, T<=T_floor

      }
    }
  } // End of loop over Meshblock

  return;
}

/*------------------------------- Frame shift ---------------------------*/
void frame_shift(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  //* To make sure that this runs only once
  if (time_last==time){
    return;
  }
  else{
    time_last = time;
  }

  Real local_cold_mass  = 0.0;
  Real global_cold_mass = 0.0;

  Real global_v3_sum = 0.0;
  Real local_v3_sum  = 0.0;

  Real front_posn_new = 0.0;
  Real v_shift_t_new = 0.0;

  Real sim_dt = pmb->pmy_mesh->dt;

  Real g = pmb->peos->GetGamma();
  Real c_s = sqrt( g*T_floor/(mu*KELVIN) );
  Real c_s_cap = 100.0;


  //* Calculate local meshblock cold mass
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        Real rho = prim(IDN,k,j,i);
        Real prs = prim(IPR,k,j,i);

        Real temp = (prs / rho) * KELVIN * mu ;

        // if (temp <= 1.1*T_floor){
        if (temp <= T_cold){
          local_cold_mass += rho*pmb->pcoord->GetCellVolume(k,j,i);
        }

        local_v3_sum += prim(IVZ,k,j,i);

      }
    }
  } // End of loop over Meshblock

  //* Calculate total cold mass
  MPI_Allreduce(&local_cold_mass, &global_cold_mass, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v3_sum, &global_v3_sum, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);

  Real v3_avg = global_v3_sum/Real(nx1*nx2*nx3);

  printf("\n_______________________________________\n");
  printf("global_cold_mass = %lf \n"  , global_cold_mass);
  printf("local_cold_mass  = %lf \n\n", local_cold_mass);

  printf("Meshblock: %ld \n", pmb->gid);
  printf("time     : %f \n", time);
  printf("\n_______________________________________\n");


  //* New front position
  front_posn_new = x3min + global_cold_mass/(L1*L2 * amb_rho*T_hot/T_floor);

  printf("L1             = %lf \n"  , L1);
  printf("L2             = %lf \n"  , L2);
  printf("x3min          = %lf \n\n", x3min);


  printf("front_posn_old = %lf \n"  , front_posn_old);
  printf("front_posn_new = %lf \n\n", front_posn_new);

  //* Required shift velocity
  v_shift_t_new = -2.0 * (front_posn_new-front_posn_old)/sim_dt;
  // Factor of two is because this is applied only once per timestep
  // instead of two

  printf("pmy_mesh->time = %lf   \n", pmb->pmy_mesh->time);
  printf("time           = %lf   \n", time);
  printf("dt             = %lf \n\n", sim_dt);

  printf("v_shift_t_new  = %lf \n\n", v_shift_t_new );
  printf("v_z            = %lf \n\n", v3_avg);



  //* If shift velocity is too high
  if (abs(v_shift_t_new) > (c_s_cap*c_s)){  // Comparison with cold gas sound speed
    v_shift_t_new /=  abs(v_shift_t_new); //Calculate sign of v_shift
    v_shift_t_new *=  c_s*c_s_cap;

    printf("_______________________________________\n");
    printf("v_shift > v_cap !... \n");
    printf("New v_shift: %lf \n", v_shift_t_new);
    printf("_______________________________________\n");
  }
  front_velocity = -1.0*v_shift_t_new;
  // printf("front_velocity = %lf \n\n", front_velocity);
  // printf("_______________________________________\n");

  // int il = pmb->is;
  int il = pmb->is - NGHOST;
  // int iu = pmb->ie;
  int iu = pmb->ie + NGHOST;

  // int jl = pmb->js;
  int jl = pmb->js- NGHOST;
  // int ju = pmb->je;
  int ju = pmb->je + NGHOST;

  // int kl = pmb->ks;
  int kl = pmb->ks- NGHOST;
  // int ku = pmb->ke;
  int ku = pmb->ke + NGHOST;

  //* Add the shift velocity
  if (time <= shift_start){
    front_posn_old = front_posn_new;
  }
  else {
    printf("\n Start of shifting... \n");
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {


          //! No issue occurs if this line is removed
          // v_shift_t_new = 0.0 ; // c_s*c_s_cap ;
          cons(IM3,k,j,i) += prim(IDN,k,j,i) * v_shift_t_new;

          Real dE_kin = SQR(prim(IVZ,k,j,i)+v_shift_t_new) - SQR(prim(IVZ,k,j,i));

          // Add the residual energy to total energy
          cons(IEN,k,j,i) += 0.5 * prim(IDN,k,j,i) * dE_kin ;

        }
      }
    } // End of loop over Meshblock
  }

  // Update front position
  // front_posn_old = front_posn_new;

  // Update shift velocity
  // v_shift_t_old = v_shift_t_new;

  return;
}
Real hst_total_cooling(MeshBlock *pmb, int iout) {
  if(pmb->lid == 0)
    return total_cooling;
  else
    return 0;
}

Real hst_front_velocity(MeshBlock *pmb, int iout) {
  if(pmb->lid == 0)
    return abs(front_velocity);
  else
    return 0.0;
}

Real hst_debug(MeshBlock *pmb, int iout) {
  if(pmb->lid == 0)
    return debug_hst;
  else
    return 0;
}

/*-------------------------------- Source ---------------------------*/
void Source(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  if (cooling_flag!=0){

    if (!cooling_flag_print_count){
      printf("___________________________________________\n");
      printf("!! Cooling included .......................\n");
      printf("___________________________________________\n");

      cooling_flag_print_count = true;
    }

    townsend_cooling(pmb, time, dt,
            prim, prim_scalar,
            bcc, cons,
            cons_scalar);

  }

  if (shift_flag!=0){

    if (!shift_flag_print_count){
      printf("___________________________________________\n");
      printf("!! Front shift included ...................\n");
      printf("___________________________________________\n");

      shift_flag_print_count = true;
    }

    frame_shift(pmb, time, dt,
            prim, prim_scalar,
            bcc, cons,
            cons_scalar);

  }

  return;

}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  vflow = pin->GetReal("problem","vflow");
  iprob = pin->GetInteger("problem","iprob");
  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem","four_pi_G");
    Real eps = pin->GetOrAddReal("problem","grav_eps", 0.0);
    SetFourPiG(four_pi_G);
    SetGravityThreshold(eps);
  }
  read_input(pin);
  cooler = make_unique<Cooling>();
  if (MAGNETIC_FIELDS_ENABLED) {

    AllocateUserHistoryOutput(12);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, Pth_sum, "Pth_sum");
    EnrollUserHistoryOutput(4, PB_sum, "PB_sum");
    EnrollUserHistoryOutput(5, Bx_sum, "Bx_sum");
    EnrollUserHistoryOutput(6, By_sum, "By_sum");
    EnrollUserHistoryOutput(7, Bz_sum, "Bz_sum");
    EnrollUserHistoryOutput(8, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(9, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(10, hst_front_velocity, "front_velocity", UserHistoryOperation::max);
    EnrollUserHistoryOutput(11, hst_debug, "total_debug");
  }
  else {

    AllocateUserHistoryOutput(7);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(4, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(5, hst_front_velocity, "front_velocity", UserHistoryOperation::max);
    EnrollUserHistoryOutput(6, hst_debug, "total_debug");

  } 
  if (adaptive) {
    threshold = pin->GetReal("problem", "thr");
    EnrollUserRefinementCondition(RefinementCondition);
  }
  if (iprob == 4 && NSCALARS > 0) {
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, PassiveDyeEntropy, "tot-S");
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholtz test

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::int64_t iseed = -1 - gid;
  Real gm1 = peos->GetGamma() - 1.0;
  read_input(pin);

  //--- iprob=1.  Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two slip-surfaces from background medium with d=1
  // moving at (+vflow), random perturbations.  This is the classic, unresolved K-H test.

  if (iprob == 1) {
    // Read problem parameters
    Real drat = pin->GetReal("problem","drat");
    Real amp = pin->GetReal("problem","amp");
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = rho_hot;
          phydro->u(IM1,k,j,i) = -(vflow + amp*(ran2(&iseed) - 0.5));
          phydro->u(IM2,k,j,i) = vflow * sin(4*PI*pcoord->x2v(i))+amp*(ran2(&iseed) - 0.5);
          phydro->u(IM3,k,j,i) = 0.0;
          if (pcoord->x2v(j) < (L2/2.0)) {
            phydro->u(IDN,k,j,i) = rho_cold;
            phydro->u(IM1,k,j,i) = vflow + amp*(ran2(&iseed) - 0.5);
            //phydro->u(IM2,k,j,i) = drat*amp*(ran2(&iseed) - 0.5);
          }
          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=2. Uniform density medium moving at +/-vflow seperated by a single shear
  // layer with tanh() profile at y=0 with a single mode perturbation, reflecting BCs at
  // top/bottom.  Based on Frank et al., ApJ 460, 777, 1996.

  if (iprob == 2) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem", "amp");
    Real a = 0.02;
    Real sigma = 0.2;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = vflow*std::tanh((pcoord->x2v(j))/a);
          phydro->u(IM2,k,j,i) = amp*std::cos(TWO_PI*pcoord->x1v(i))
                                 *std::exp(-(SQR(pcoord->x2v(j)))/SQR(sigma));
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=3.  Test in SR paper (Beckwith & Stone, ApJS 193, 6, 2011).  Gives two
  // resolved shear layers with tanh() profiles for velocity and density located at
  // y = +/- 0.5, density one in middle and 0.01 elsewhere, single mode perturbation.

  if (iprob == 3) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    Real a = 0.01;
    Real sigma = 0.1;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 0.505 + 0.495
                                 *std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM1,k,j,i) = vflow*std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM2,k,j,i) =
              amp*vflow*std::sin(TWO_PI*pcoord->x1v(i))
              *std::exp(-((std::abs(pcoord->x2v(j))-0.5)
                          *(std::abs(pcoord->x2v(j))-0.5))/(sigma*sigma));
          if (pcoord->x2v(j) < 0.0) phydro->u(IM2,k,j,i) *= -1.0;
          phydro->u(IM1,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM2,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=4.  "Lecoanet" test, resolved shear layers with tanh() profiles for velocity
  // and density located at z1=0.5, z2=1.5 two-mode perturbation for fully periodic BCs

  // To promote symmetry of FP errors about midplanes, rescale z' = z - 1. ; x' = x - 0.5
  // so that domain x1 = [-0.5, 0.5] and x2 = [-1.0, 1.0] is centered about origin
  if (iprob == 4) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    // unstratified problem is the default
    Real drho_rho0 = pin->GetOrAddReal("problem", "drho_rho0", 0.0);
    // set background vx to nonzero to evolve the KHI in a moving frame
    Real vboost = pin->GetOrAddReal("problem", "vboost", 0.0);
    Real P0 = 10.0;
    Real a = 0.05;
    Real sigma = 0.2;
    // Initial condition's reflect-and-shift symmetry, x1-> x1 + 1/2, x2-> -x2
    // is preserved in new coordinates; hence, the same flow is solved twice in this prob.
    Real z1 = -0.5;  // z1' = z1 - 1.0
    Real z2 = 0.5;   // z2' = z2 - 1.0

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          // Lecoanet (2015) equation 8a)
          Real dens = 1.0 + 0.5*drho_rho0*(std::tanh((pcoord->x2v(j) - z1)/a) -
                                           std::tanh((pcoord->x2v(j) - z2)/a));
          phydro->u(IDN,k,j,i) = dens;

          Real v1 = vflow*(std::tanh((pcoord->x2v(j) - z1)/a)
                           - std::tanh((pcoord->x2v(j) - z2)/a) - 1.0) // 8b)
                    + vboost;
          // Currently, the midpoint approx. is applied in the momenta and energy calc
          phydro->u(IM1,k,j,i) = v1*dens;

          // NOTE ON FLOATING-POINT SHIFT SYMMETRY IN X1:
          // There is no scaling + translation coordinate transformation that would
          // naturally preserve this symmetry when calculating x1 coordinates in
          // floating-point representation. Need to correct for the asymmetry of FP error
          // by modifying the operands.  Centering the domain on x1=0.0 ensures reflective
          // symmetry, x1' -> -x1 NOT shift symmetry, x1' -> x1 + 0.5 (harder guarantee)

          // For example, consider a cell in the right half of the domain with x1v > 0.0,
          // so that shift symmetry should hold w/ another cell's x1v'= -0.5 + x1v < 0.0

          // ISSUE: sin(2*pi*(-0.5+x1v)) != -sin(2*pi*x1v) in floating-point calculations
          // The primary FP issues are twofold: 1) different rounding errors in x1v, x1v'
          // and 2) IEEE-754 merely "recommends" that sin(), cos(), etc. functions are
          // correctly rounded. Note, glibc library doesn't provide correctly-rounded fns

          // 1) Since x1min = -0.5 can be perfectly represented in binary as -2^{-1}:
          // double(x1v')= double(double(-0.5) + double(x1v)) = double(-0.5 + double(x1v))
          // Even if x1v is also a dyadic rational -> has exact finite FP representation:
          // x1v'= double(-0.5 + double(x1v)) = double(-0.5 + x1v) ?= (-0.5 + x1v) exactly

          // Sterbenz's Lemma does not hold for any nx1>4, so cannot guarantee exactness.
          // However, for most nx1 = power of two, the calculation of ALL cell center
          // positions x1v will be exact. For nx1 != 2^n, differences are easily observed.

          // 2) Even if the rounding error of x1v (and hence x1v') is zero, the exact
          // periodicity of trigonometric functions (even after range reduction of input
          // to [-pi/4, pi/4], e.g.) is NOT guaranteed:
          // sin(2*pi*(-0.5+x1v)) = sin(-pi + 2*pi*x1v) != -sin(2*pi*x1v)

          // WORKAROUND: Average inexact sin() with -sin() sample on opposite x1-half of
          // domain The assumption of periodic domain with x1min=-0.5 and x1max=0.5 is
          // hardcoded here (v2 is the only quantity in the IC with x1 dependence)

          Real ave_sine = std::sin(TWO_PI*pcoord->x1v(i));
          if (pcoord->x1v(i) > 0.0) {
            ave_sine -= std::sin(TWO_PI*(-0.5 + pcoord->x1v(i)));
          } else {
            ave_sine -= std::sin(TWO_PI*(0.5 + pcoord->x1v(i)));
          }
          ave_sine /= 2.0;

          // translated x1= x - 1/2 relative to Lecoanet (2015) shifts sine function by pi
          // (half-period) and introduces U_z sign change:
          Real v2 = -amp*ave_sine
                    *(std::exp(-(SQR(pcoord->x2v(j) - z1))/(sigma*sigma)) +
                      std::exp(-(SQR(pcoord->x2v(j) - z2))/(sigma*sigma))); // 8c), mod.
          phydro->u(IM2,k,j,i) = v2*dens;

          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                                 + SQR(phydro->u(IM2,k,j,i))
                                                 + SQR(phydro->u(IM3,k,j,i)) )
                                   /phydro->u(IDN,k,j,i);
          }
          // color concentration of passive scalar
          if (NSCALARS > 0) {
            Real concentration = 0.5*(std::tanh((pcoord->x2v(j) - z2)/a)  // 8e)
                                      - std::tanh((pcoord->x2v(j) - z1)/a) + 2.0);
            // uniformly fill all scalar species to have equal concentration
            constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;
            for (int n=0; n<NSCALARS; ++n) {
              pscalars->s(n,k,j,i) = 1.0/scalar_norm*concentration*phydro->u(IDN,k,j,i);
            }
          }
        }
      }
    }
    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem", "b0");
      b0 = b0/std::sqrt(4.0*(PI));
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0*std::tanh((std::abs(pcoord->x2v(j)) - 0.5)/a);
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*SQR(pfield->b.x1f(k,j,i));
            }
          }
        }
      }
    }
  }

  //--- iprob=5. Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two resolved slip-surfaces from background medium
  // with d=1 moving at (+vflow), with m=2 perturbation, for the AMR test.

  if (iprob == 5) {
    // Read problem parameters
    Real a = pin->GetReal("problem","a");
    Real sigma = pin->GetReal("problem","sigma");
    Real drat = pin->GetReal("problem","drat");
    Real amp = pin->GetReal("problem","amp");
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real w=(std::tanh((std::abs(pcoord->x2v(j))-0.25)/a)+1.0)*0.5;
          phydro->u(IDN,k,j,i) = w+(1.0-w)*drat;
          phydro->u(IM1,k,j,i) = w*vflow-(1.0-w)*vflow*drat;
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*amp
                                 * std::sin(2.0*TWO_PI*pcoord->x1v(i))
                                 * std::exp(-SQR(std::abs(pcoord->x2v(j))-0.25)
                                            /(sigma*sigma));
          phydro->u(IM3,k,j,i) = 0.0;
          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.25*(SQR(phydro->u(IM1,k,j,i)) +
                                SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  return;
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin){
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){

  Real g = peos->GetGamma();

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {

        Real lum_cell = 0.0;
        Real temp = (phydro->u(IPR,k,j,i) / phydro->u(IDN,k,j,i)) * KELVIN * mu ;

        if (temp > T_floor) {

            //* For Max's townsend cooling

            if (temp<T_cut){
              Real rho      = phydro->u(IDN,k,j,i);
              Real temp_cgs = temp;
              Real rho_cgs  = rho * unit_density;
              Real dt_cgs   = pmy_mesh->dt * unit_time;
              Real cLfac    = Lambda_fac;

              Real temp_new = fmax(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, cLfac), T_floor);

              Real ccool_2 = ((temp_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);

              lum_cell-= ccool_2;
            }

            //*_____________________________

        }

        else{  // If T<=T_floor

          Real ccool_2 = ((T_floor-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);
          // lum_cell-= ccool_2;

        }

        user_out_var(0,k,j,i) = lum_cell/pmy_mesh->dt;

      }
    }
  }
}
// refinement condition: velocity gradient

int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vgmax = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real vgy = std::abs(w(IVY,k,j,i+1) - w(IVY,k,j,i-1))*0.5;
        Real vgx = std::abs(w(IVX,k,j+1,i) - w(IVX,k,j-1,i))*0.5;
        Real vg  = std::sqrt(vgx*vgx + vgy*vgy);
        if (vg > vgmax) vgmax = vg;
      }
    }
  }
  if (vgmax > threshold) return 1;
  if (vgmax < 0.5*threshold) return -1;
  return 0;
}

namespace {
Real PassiveDyeEntropy(MeshBlock *pmb, int iout) {
  Real total_entropy = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &r = pmb->pscalars->r;
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=is; i<=ie; i++) {
        // no loop over NSCALARS; hardcode assumption that NSCALARS=1
        Real specific_entropy = -r(0,k,j,i)*std::log(r(0,k,j,i));
        total_entropy += volume(i)*w(IDN,k,j,i)*specific_entropy;  // Lecoanet (2016) eq 5
      }
    }
  }
  return total_entropy;
}
} // namespace
