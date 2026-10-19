#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// Simple convection climate model
// Conventions: Used Units are SI Units
// Atmosphere lists

// Definition of constants
#define G 9.80665
#define CP 1004.0
#define R_DRY 287.05 //specific gas constant for dry air (J/kg /K)
#define KAPPA (2.0/7.0)
#define GAMMA_WET 0.0065 // mean atmospheric temperature gradient (K/m)
#define GAMMA_DRY 0.01 // dry adiabatic temperature gradient (K/m)  
#define SIGMA 5.670374419e-8 // Stefan-Boltzmann constant (W/m^2 K^4)


// Definition of useful helper functions
double TtoTheta (double p, double T){ // an inverse calc could be useful and implemented here too
    /*Convert a given temperature to the respective potential temperature or vice-versa, p in Pa, T in K
    For calculation o f T from theta, the T argument has to hold the theta value  */
    double p0 = 100000;

    double theta = T * pow(p0/p, KAPPA);
    double temp = T * pow(p/p0, KAPPA); // not yet implemented: simple if statement on additional arg

    // for calculation of T from theta: assume T is the theta value
    return theta;
}

void init_heights_temps(const double p[], double T[], double z[],
                        double gamma, int num_levels) {
    /* Initializes temperature T[] (K) and height z[] (m) on pressure levels p[] (Pa).
       Ordering: p[0] = 0 Pa (top), p[num_levels-1] = ground pressure.
       T[num_levels-1] must be set by the caller as the surface temperature.
       Heights are relative to the ground (z = 0 at the last level). */

    const int s = num_levels - 1;            // surface index
    const double T0 = T[s];
    const double p0 = p[s];
    const double exponent = R_DRY * gamma / G;

    z[s] = 0.0;

    for (int k = 0; k < s; k++) {// loop over all levels above the ground
        T[k] = T0 * pow(p[k] / p0, exponent);
        z[k] = (T0 - T[k]) / gamma;
    }
}

double heating_rate(double Eabs, double dt, double dp) {
    return Eabs * dt * G / (dp * CP);
}

int main(int argc, char **argv)
{
    // Model constraints
    int nlayer = 20;
    int n_timesteps = 1000;

    // Define atmosphere
    double Eabs = 240;
    double dt = 3600; // time increment for model (s) TODO: check unit
    double p_ref = 100000.0; // reference pressure for potential temperature calculation
    double T_lowest_init = 288; // initial temperature of lowest layer (K)

    double p_lev[nlayer + 1]; // pressure levels (Pa)
    double theta_lay[nlayer];
    double T_lay[nlayer]; // layer temperatures (K)
    double T_lev[nlayer + 1];   // level temperatures (K)

    double dp = 100000.0 / nlayer; // pressure difference betw. levels (Pa)
    double p_lay[nlayer]; // mean assumed pressure for each layer (Pa)
    double z_lev[nlayer + 1]; // geometric height levels (m)
    double z_lay[nlayer]; // geometric height of each layer (m)

    double E_up [nlayer+1]; // upward energy flux at each level (W/m^2)
    double E_down [nlayer+1]; // downward energy flux at each level (W/m^2)

    // Initialize pressure levels
    for (int i = 0; i < nlayer + 1; i++) //
    {
        p_lev[i] = i * dp; // list of level pressures
    }

    for (int i = 0; i < nlayer; i++) //initialzie layer pressures
    {
        p_lay[i] = (p_lev[i] + p_lev[i+1]) / 2.0; // pressure for each layer as mean of the two pressure levels
    }

    // initialize temperature profile and geometric heights
    T_lev[nlayer] = T_lowest_init; // set initial temperature on the lowest level
    init_heights_temps(p_lev, T_lev, z_lev, GAMMA_DRY, nlayer+1); // z used as helper for temp init for now

    for (int i = 0; i < nlayer; i++) 
    {
        T_lay[i] = (T_lev[i] + T_lev[i+1]) / 2.0; // layer temperatures as mean of the two level temperatures
        z_lay[i] = (z_lev[i] + z_lev[i+1]) / 2.0; // geometric height of each layer as mean of the two level heights
    }


    // calculate potential temperature
    for (int i = 0; i < nlayer; i++)
    {   
        theta_lay[i] = TtoTheta(p_lay[i], T_lay[i]); // calculate potential temperature for each layer given the starting temp an layer press
    }


    // Calculate temperature change from energy input
    double dT = heating_rate(Eabs, dt, dp); // calculate temperature change from energy input 

    // Print results
    for (int i = 0; i < nlayer; i++)
    {
        printf("Layer %d: p = %.1f Pa, theta = %.2f K, temperature = %.2f K, altitude = %.2f m\n ",
               i, p_lay[i], theta_lay[i], T_lay[i], z_lay[i]);
    }

    printf("Temperature change: %.6f K\n", dT);

    int count = 0; 
    while (count < n_timesteps) {
        count++;
        
        //T[nlayer-1] += dT; // Update temperature of the lowermost layer
        //theta[nlayer-1] = T[nlayer-1] * pow(p_ref / play[nlayer-1], 0.286);
        // sort by theta and update T accordingly
        // radiative transfer from T 


    return 0;
}