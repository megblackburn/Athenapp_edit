#ifndef COOLING_CLASS
#define COOLING_CLASS

#include "../athena.hpp"
#include "../athena_arrays.hpp"

static bool DEBUG_FLAG = false;

class Cooling {
public:
  Cooling();
  ~Cooling();

  Real tfloor;
  Real townsend(Real temp, Real rho, Real const dt, Real Lambda_fac = 1.0);
  Real Lambda(Real temp);
  Real tcool(Real temp, Real rho);

private:
  Real const mh = 1.6605e-24;
  Real const kb = 1.380648e-16;

  int nbins;
  
  Real const_factor;

  AthenaArray<Real> cool_t;
  AthenaArray<Real> cool_tef;
  AthenaArray<Real> cool_coef;
  AthenaArray<Real> cool_index;
};

Cooling::Cooling() {

  if (DEBUG_FLAG){
    printf("Start of Cooling::Cooling!\n");
  }

  // Compute the constant factor needed to compute new TEFs
  Real g  = 5.0/3.0; // adiabatic index
  Real X = 0.7; auto Z = 0.02; // Fully ionized, solar abundances
  Real mu   = 1.0/(2.0*X + 0.75*(1.0-X-Z) + Z/2.0);
  Real mu_e = 2.0/(1.0+X);
  Real mu_h = 1.0/X;
  const_factor = (1.0e-23)*(g-1.0)*mu/(kb*mu_e*mu_h*mh);

  // Initialize the cooling table
  nbins = 40;

  cool_t.NewAthenaArray(nbins);
  cool_tef.NewAthenaArray(nbins);
  cool_coef.NewAthenaArray(nbins);
  cool_index.NewAthenaArray(nbins);

  // Temperatures in [K]
  cool_t(0)  = 10000.0;
  cool_t(1)  = 12600.0;
  cool_t(2)  = 15890.0;
  cool_t(3)  = 20020.0;
  cool_t(4)  = 25240.0;
  cool_t(5)  = 31810.0;
  cool_t(6)  = 40090.0;
  cool_t(7)  = 50530.0;
  cool_t(8)  = 63680.0;
  cool_t(9)  = 80260.0;
  cool_t(10) = 101200.0;
  cool_t(11) = 127500.0;
  cool_t(12) = 160700.0;
  cool_t(13) = 202600.0;
  cool_t(14) = 255300.0;
  cool_t(15) = 321800.0;
  cool_t(16) = 405600.0;
  cool_t(17) = 511100.0;
  cool_t(18) = 644200.0;
  cool_t(19) = 812000.0;
  cool_t(20) = 1000000.0;
  cool_t(21) = 1259000.0;
  cool_t(22) = 1585000.0;
  cool_t(23) = 1995000.0;
  cool_t(24) = 2512000.0;
  cool_t(25) = 3162000.0;
  cool_t(26) = 3981000.0;
  cool_t(27) = 5012000.0;
  cool_t(28) = 6310000.0;
  cool_t(29) = 7943000.0;
  cool_t(30) = 10000000.0;
  cool_t(31) = 12590000.0;
  cool_t(32) = 15850000.0;
  cool_t(33) = 19950000.0;
  cool_t(34) = 25120000.0;
  cool_t(35) = 31620000.0;
  cool_t(36) = 39810000.0;
  cool_t(37) = 50120000.0;
  cool_t(38) = 63100000.0;
  cool_t(39) = 79430000.0;

  // Cooling Coefficient [1e-23 ergs*cm3/s]
  cool_coef(0) = 1.6408984689285624;
  cool_coef(1) = 5.789646575948292;
  cool_coef(2) = 18.797203755396648;
  cool_coef(3) = 16.7384754689852;
  cool_coef(4) = 11.274384717759935;
  cool_coef(5) = 9.95038422958871;
  cool_coef(6) = 11.302144847043829;
  cool_coef(7) = 15.819149070534786;
  cool_coef(8) = 25.224636283348048;
  cool_coef(9) = 38.02107089248533;
  cool_coef(10) = 43.98219098299675;
  cool_coef(11) = 41.277704007796586;
  cool_coef(12) = 41.95311185975414;
  cool_coef(13) = 45.260670345801;
  cool_coef(14) = 47.275626188961176;
  cool_coef(15) = 32.21420131907784;
  cool_coef(16) = 24.350976818250636;
  cool_coef(17) = 23.383616480583676;
  cool_coef(18) = 18.333394532081098;
  cool_coef(19) = 14.89691888284402;
  cool_coef(20) = 14.392505898454834;
  cool_coef(21) = 13.027915287005817;
  cool_coef(22) = 11.671262753284271;
  cool_coef(23) = 9.070904785425046;
  cool_coef(24) = 6.489695397654223;
  cool_coef(25) = 4.766239129792971;
  cool_coef(26) = 3.7811870710765074;
  cool_coef(27) = 3.313622783657129;
  cool_coef(28) = 3.0600313080475674;
  cool_coef(29) = 2.9993768457216112;
  cool_coef(30) = 2.9491332141250552;
  cool_coef(31) = 2.744653611808266;
  cool_coef(32) = 2.3449511265716;
  cool_coef(33) = 2.0169621177549892;
  cool_coef(34) = 1.8907205849384978;
  cool_coef(35) = 1.91584885606706;
  cool_coef(36) = 2.056870288868004;
  cool_coef(37) = 2.233680315878366;
  cool_coef(38) = 2.4097186710383474;
  cool_coef(39) = 2.5502102007949023;

  // Cooling power index
  cool_index(0) = 5.455488390256632;
  cool_index(1) = 5.076170519863754;
  cool_index(2) = -0.5020655826640291;
  cool_index(3) = -1.7055659800651979;
  cool_index(4) = -0.5399688186820728;
  cool_index(5) = 0.550609170202909;
  cool_index(6) = 1.4527662908446985;
  cool_index(7) = 2.0172644735605223;
  cool_index(8) = 1.773197476674277;
  cool_index(9) = 0.6282445620956022;
  cool_index(10) = -0.2747076405016009;
  cool_index(11) = 0.07013182420220869;
  cool_index(12) = 0.32752568568776125;
  cool_index(13) = 0.1883881016798681;
  cool_index(14) = -1.6570303755459093;
  cool_index(15) = -1.209120245966656;
  cool_index(16) = -0.17533183860418153;
  cool_index(17) = -1.0512755674245657;
  cool_index(18) = -0.896664392554265;
  cool_index(19) = -0.16540667885641686;
  cool_index(20) = -0.43250361812273735;
  cool_index(21) = -0.4775539072045259;
  cool_index(22) = -1.0956186284443203;
  cool_index(23) = -1.453147878451421;
  cool_index(24) = -1.3412596915753237;
  cool_index(25) = -1.0051719479026813;
  cool_index(26) = -0.573142729390977;
  cool_index(27) = -0.3457087236213044;
  cool_index(28) = -0.08698732111048613;
  cool_index(29) = -0.07335511773234596;
  cool_index(30) = -0.3119882060952377;
  cool_index(31) = -0.6835132944311395;
  cool_index(32) = -0.6549261784681947;
  cool_index(33) = -0.2804886559029823;
  cool_index(34) = 0.05737205818565948;
  cool_index(35) = 0.30836313806582183;
  cool_index(36) = 0.3580735000106496;
  cool_index(37) = 0.3293929876114671;
  cool_index(38) = 0.24620665148692336;
  cool_index(39) = 0.10953503955831644;

  // Calculate All TEFs Y_k recursively (Eq. A6)
  cool_tef(nbins-1) = 0.0; // Last Y_N = 0
  for (int i=nbins-2; i>=0; i--) {
    Real t_n    = cool_t(nbins-1);
    Real coef_n = cool_coef(nbins-1);

    Real t_i   = cool_t(i);
    Real tef_i = cool_tef(i);
    Real coef  = cool_coef(i);
    Real slope = cool_index(i);

    Real sm1  = slope - 1.0;
    Real step = (coef_n/coef)*(t_i/t_n)*(std::pow(t_i/cool_t(i+1),sm1)-1)/sm1;
    cool_tef(i) = cool_tef(i+1) - step;
  }

  // Set temperature floor
  tfloor = cool_t(0);


  if (DEBUG_FLAG){
    printf("End of Cooling::Cooling!\n");
  }

}

Cooling::~Cooling() {
  cool_t.DeleteAthenaArray();
  cool_tef.DeleteAthenaArray();
  cool_coef.DeleteAthenaArray();
  cool_index.DeleteAthenaArray();
}


Real Cooling::townsend(Real temp, Real rho, Real const dt, Real Lambda_fac)
{
if (DEBUG_FLAG){
    printf("Start of Cooling::townsend!\n");
  }

  if (DEBUG_FLAG){
    printf("Check that temperature is above the cooling floor\n");
  }
  // Check that temperature is above the cooling floor
  // if (temp < tfloor) return tfloor;


  if (DEBUG_FLAG){
    printf("Get Reference values from the last bin\n");
  }
  // Get Reference values from the last bin
  Real t_n    = cool_t(nbins-1);
  Real coef_n = cool_coef(nbins-1) * Lambda_fac;

  if (DEBUG_FLAG){
    printf("Get the index of the right temperature bin\n");
  }
  // Get the index of the right temperature bin
  int idx = 0;
  while ((idx < nbins-2) && (cool_t(idx+1) < temp)) { idx += 1; }

  if (DEBUG_FLAG){
    printf("Look up the corresponding slope and coefficient\n");
  }
  // Look up the corresponding slope and coefficient
  Real t_i   = cool_t(idx);
  Real tef_i = cool_tef(idx);
  Real coef  = cool_coef(idx) * Lambda_fac;
  Real slope = cool_index(idx);

  if (DEBUG_FLAG){
    printf("Compute the Temporal Evolution Function Y(T) (Eq. A5)\n");
  }
  // Compute the Temporal Evolution Function Y(T) (Eq. A5)
  Real sm1 = slope - 1.0;
  Real tef = tef_i + (coef_n/coef)*(t_i/t_n)*(std::pow(t_i/temp,sm1)-1)/sm1;

  if (DEBUG_FLAG){
    printf("Compute the adjusted TEF for new timestep (Eqn. 26)\n");
  }
  // Compute the adjusted TEF for new timestep (Eqn. 26)
  Real tef_adj = tef + rho*coef_n*const_factor*dt/t_n;

  if (DEBUG_FLAG){
    printf("TEF is a strictly decreasing function and new_tef > tef\n");
  }
  // TEF is a strictly decreasing function and new_tef > tef
  // Check if the new TEF falls into a lower bin
  // If so, update slopes and coefficients
  while ((idx > 0) && (tef_adj > cool_tef(idx))) {
    idx -= 1;
    t_i   = cool_t(idx);
    tef_i = cool_tef(idx);
    coef  = cool_coef(idx) * Lambda_fac;  //* Later correction
    slope = cool_index(idx);
  }

  if (DEBUG_FLAG){
    printf("Compute the Inverse Temporal Evolution Function Y^{-1}(Y) (Eq. A7)\n");
  }
  // Compute the Inverse Temporal Evolution Function Y^{-1}(Y) (Eq. A7)
  Real oms  = 1.0 - slope;
  Real tnew = t_i*std::pow(1-oms*(coef/coef_n)*(t_n/t_i)*(tef_adj-tef_i),1/oms);

  if (DEBUG_FLAG){
    printf("End of Cooling::townsend!\n");
  }
  // Return the new temperature if it is still above the temperature floor
  //return std::max(tnew,tfloor);
  return tnew;
}


Real Cooling::Lambda(Real temp)
{
  if (DEBUG_FLAG){
    printf("Start of Cooling::Lambda!\n");
  }

  if(temp<cool_t(0) ||(temp> cool_t(nbins-1)))
    return NAN;

  int idx = 0;
  while ((idx<nbins-2) && (cool_t(idx+1)<temp)) {idx += 1; }

  Real t_i = cool_t(idx);
  Real coef = cool_coef(idx);
  Real slope = cool_index(idx);

  if (DEBUG_FLAG){
    printf("End of Cooling::Lambda!\n");
  }
  return coef * pow(temp/t_i, slope);
}

Real Cooling::tcool(Real temp, Real rho)
{
  if (DEBUG_FLAG){
    printf("Start of Cooling::tcool!\n");
  }
  if((temp<cool_t(0)) || (temp > cool_t(nbins-1)))
    return NAN;
  Real L = Cooling::Lambda(temp);

  if (DEBUG_FLAG){
    printf("End of Cooling::tcool!\n");
  }
  return temp / (const_factor * rho * L);

}

#endif
