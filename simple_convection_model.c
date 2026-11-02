#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// Simple convection climate model
// Conventions: Used Units are SI Units
// Atmosphere lists are ordered from 0 (top) to nlayer/nlayer-1 for levels and layers:

//  0-------------------0
//            0
//  1-------------------1
//            1
//
//               ...
//  nlay-1-------------------nlay-1
//              nlay-2
//    nlay-------------------nlay
//              nlay-1


  ////////////////
 /// TODO List //
////////////////

// O Add dynamic stopping condition when in equilibrium (Check temp change over time < epsilon)
// O Clear the alpha/tau confusion
// O Implement Radiances as in Lecture



  ////////////////////////////// 
 /// Definition of constants //
//////////////////////////////

#define G 9.80665
#define CP 1004.0
#define R_DRY 287.05 //specific gas constant for dry air (J/kg /K)
#define KAPPA (R_DRY/CP) // exponent of potential temperature calculation
#define GAMMA_WET 0.0065 // mean atmospheric temperature gradient (K/m)
#define GAMMA_DRY 0.01 // dry adiabatic temperature gradient (K/m)  
#define SIGMA 5.670374419e-8 // Stefan-Boltzmann constant (W/m^2 K^4)
#define P_REF 100000.0 // reference pressure for potential temperature calculation
#define PI 3.14159265358979323846 // ensure consistent PI - could likely use M_PI too

  ////////////////////////////////////////////
 /// Definition of useful helper functions //
////////////////////////////////////////////

double TtoTheta (double p, double T, bool inverse){ // an inverse calc could be useful and implemented here too
    /*Convert a given temperature to the respective potential temperature or vice-versa, p in Pa, T in K
    For calculation o f T from theta, the T argument has to hold the theta value  */

    double theta = T * pow(P_REF/p, KAPPA);
    double temp = T * pow(p/P_REF, KAPPA); 

    // for calculation of T from theta: assume T is the theta value
    if (inverse) {
        return temp;
    } else {
        return theta;
    }
}

void init_heights_temps(double p[], double T[], double z[],
                        double gamma, int num_levels) {
    /* Initializes temperature T[] (K) and height z[] (m) on pressure levels p[] (Pa).
       Ordering: p[0] = 0 Pa (top), p[num_levels-1] = ground pressure.
       T[num_levels-1] must be set by the caller as the surface temperature.
       Heights are relative to the ground (z = 0 at the last level). */

    int s = num_levels - 1;            // surface index
    double T0 = T[s];
    double p0 = p[s];
    double exponent = R_DRY * gamma / G;

    z[s] = 0.0;

    for (int k = 0; k < s; k++) {// loop over all levels above the ground, hypsometric eq.
        T[k] = T0 * pow(p[k] / p0, exponent);
        z[k] = (T0 - T[k]) / gamma;
    }
}

double heating_rate(double Eabs, double dt, double dp) {
    return Eabs * dt * G / (dp * CP);
}

void convection(double T_lay[], double theta_lay[], double p_lay[], int nlayer) {

    int unsorted = 1; // flag to indicate that the layers are unsorted
    while (unsorted == 1){
        unsorted = 0;
        for (int i = nlayer-1; i > -1 ; i--){ // sort from bottom to top for faster runtime
            if (i == 0){} // this case is covered by i = 1 already
            else if (theta_lay[i-1] < theta_lay[i]){ // check for convective instability 
                unsorted = 1;//code reaches this area!
                // swap potential temperatures of the two layers

                double temp_theta = theta_lay[i]; // intermediate step for sorting 
                theta_lay[i] = theta_lay[i-1];
                theta_lay[i-1] = temp_theta;
                }
            }
        }
    for (int i = 0; i < nlayer; i++) {
        T_lay[i] = TtoTheta(p_lay[i], theta_lay[i], 1);
    }
}

// calculates the spectral radiance for given temperature and wavelength
double planck(double T, double lambda){// TODO: potentially already sum over range of lambdas here? Alternatively: Loop calling the func for mult. lambda
    double h = 6.62607015e-34;         // TODO: check lambda units! Use nm and convert to m in func?
    double c = 299792458;
    double k_B = 1.380649e-23;

    double term1 = (2*h*c*c)/pow(lambda, 5);
    double term2 = 1/(exp((h*c)/(lambda*k_B*T))-1);

    return term1 * term2;
}

// calculates radiative transfer using irradiance. Modifies E_up and E_down to be used for temperature calulation
void irradiance_transfer(double T_lay[], double alphas[], double transmiss[], 
                        double E_up[], double E_down[], double Eground, int nlayer) {

    // Calculate Flux boundary conditions
    E_up[nlayer] = Eground; // flux from lowest layer is just ground flux (+ emission of layer itself which is calculated in the upward flux loop)//SIGMA * pow(T_lay[nlayer-1], 4); // lowest layer heating the layer above (Boltzmann)
    E_down[0] = 0; // highest layer does not receive significant amounts of energy from space

    //calculate downwards fluxes at each level
    for (int i = 0; i < nlayer; i++){ // so, anstatt i < nlayer +1, weil wir vom oberstem level bis zum boden E_down berechnen
        E_down[i+1] =  E_down[i] * alphas[i] + SIGMA* pow(T_lay[i], 4) * transmiss[i];
    }

    E_up[nlayer] = Eground + E_down[nlayer]; // surface balance, what is radiated into ground must be reemitted!
    // calculate upwards fluxes at each level
    for (int i = nlayer-1; i > -1 ; i--){ // TODO: check correct loop boundaries and indicecs
        E_up[i] = E_up[i+1] * alphas[i] + SIGMA * pow(T_lay[i], 4) * transmiss[i];
    }
}

// Radiative transfer via radiance (grey atmosphere, Lambertian surface).
// Integrates L over mu_j to obtain E_up and E_down (E = 2*pi * sum L(mu_j) * mu_j * dmu)
void radiance_transfer(double T_lay[], double tau[], double mu[], int nmu, double dmu,
                       double E_up[], double E_down[], double Eground, int nlayer) {

    double L[nlayer + 1];   // scratch array, overwritten for each mu_j

    // zero arrays as thei are randomly initialized (problem for += summation)
    for (int i = 0; i < nlayer + 1; i++) {
        E_up[i] = 0.0;
        E_down[i] = 0.0;
    }

    //downward
    for (int j = 0; j < nmu; j++) {
        L[0] = 0.0;// boundary condition at TOA
        for (int i = 0; i < nlayer; i++) {
            double x = tau[i] / mu[j];
            double t = exp(-x);// transmissivity for this path
            double a = -expm1(-x);
            double B = SIGMA * pow(T_lay[i], 4) / PI; // grey Planck radiance TODO: use wavelength as parameter and loop outside this func to integrate 

            L[i+1] = L[i] * t + B * a;
            E_down[i+1] += 2.0 * PI * L[i+1] * mu[j] * dmu;
        }
    }

    //surface: same radiance in all directions
    E_up[nlayer] = Eground + E_down[nlayer];
    double L_surf = E_up[nlayer] / PI;

    //upward
    for (int j = 0; j < nmu; j++) {
        L[nlayer] = L_surf;
        for (int i = nlayer - 1; i >= 0; i--) {
            double x = tau[i] / mu[j];
            double t = exp(-x);
            double a = -expm1(-x);
            double B = SIGMA * pow(T_lay[i], 4) / PI;

            L[i] = L[i+1] * t + B * a;
            E_up[i] += 2.0 * PI * L[i] * mu[j] * dmu;
        }
    }
}

// Calculates layer T and theta using E_down and E_up from the irradiance/radiance transfer function 
void adjust_temps(double E_down[], double E_up[],  double E_net[], double E_abs[], double T_lay[], double theta_lay[], 
                 double p_lay[], double dt, double dp, int nlayer){

    //calculate E_net for each level (E_net = E_dn - E_up)
    for (int i = 0; i < nlayer + 1;  i++){
        E_net[i] = E_down[i] - E_up[i]; // upwards flux is "leaving" earth -> negative sign
    }

    //calculate E_abs for each layer (E_abs = Enet von oben - E_net von unten): problematischster Teil
    for (int i = 0; i < nlayer; i++){   // Layer quantity
        E_abs[i] = E_net[i] - E_net[i+1]; // minus makes sense as fluxes pointing upwards are negative --> - (a negative quantity) > 0, basically adding energy
    }

    // update temperature for each layer using E_abs 
    for (int i = 0; i < nlayer; i++){
        double dT = heating_rate(E_abs[i], dt, dp);
        T_lay[i] += dT;
    }
    
    // update potential temperatures for each layer
    for (int i = 0; i < nlayer; i++){
        theta_lay[i] = TtoTheta(p_lay[i], T_lay[i], 0);
    }
}




  ////////////////////////////
 /// Main program function //
////////////////////////////


int main(int argc, char **argv)
{
    // Model constraints
    int nlayer = 40;
    int n_timesteps = 8000;

    // Constants defining atmosphere
    double Eground = 240; // default for incoming energy on ground (W/m^2)
    double dt = 3600; // time increment for model (s) 
    double dp = 100000.0 / nlayer; // pressure difference betw. levels (Pa)
    double T_lowest_init = 288; // iniial temperature of lowest layer (K)

    // constants defining radiation behaviour
    double alpha = 0.944; // absorption coefficient (const for all layers currently); 0.944 yield Eq. Temp of 288.55
    double transmissivity = 1-alpha;

    //atmospheric quantities
    double p_lev[nlayer + 1]; // pressure levels (Pa)
    double T_lev[nlayer + 1];   // level temperatures (K)
    double z_lev[nlayer + 1]; // geometric height levels (m)

    double theta_lay[nlayer]; // layer potential temperatures
    double T_lay[nlayer]; // layer temperatures (K)
    double p_lay[nlayer]; // mean assumed pressure for each layer (Pa)
    double z_lay[nlayer]; // geometric height of each layer (m)


    // Energy related quantities
    double E_up [nlayer+1]; // upward energy flux at each level (W/m^2) (Irradiance)
    double E_down [nlayer+1]; // downward energy flux at each level (W/m^2) (Irradiance)
    double E_net [nlayer+1]; // net energy flux at each level(W/m^2) (Irradiance)

    double L_up[nlayer+1]; // upward energy flux at each level (W/m^2 ster) (Radiance)
    double L_down[nlayer+1]; // downward energy flux at each level (W/m^2 ster) (Radiance)
    double L_net[nlayer+1]; // net energy flux at each level (W/m^2 ster) (Radiance)
    
    double E_abs [nlayer]; // LAYER QTY!! Net absorbed energy at each layer (E_net[i] - E_net[i+1]) (Irradiance)

    double L_abs [nlayer]; // LAYER QTY!! Net absorbed energy at each layer (L_net[i] - L_net[i+1]) (Radiance)

    double alphas [nlayer]; // absorption coefficients of each layer
    double transmiss [nlayer]; // transmissivity of each layer (t = 1-alpha)

    int nmu = 20; // amount of mu angles to integrate over for radiance transfer, will/can/should be increased
    double dmu = 1.0 / nmu; // step size for mu integration
    double mu[nmu]; // mu angles for radiance integration (cosine of zenith angle) 
    for (int j = 0; j < nmu; j++) mu[j] = (j + 0.5) * dmu;   // mid-point rule for integration

    double tau_total = 1.0; // total optical depth of the atmosphere (constant for now, profile TBD)
    double tau[nlayer]; // optical depth of each layer (tau = tau_total * dp / 100000.0 for now, profile TBD)
    for (int i = 0; i < nlayer; i++) tau[i] = tau_total * dp / 100000.0; // optical depth of each layer (constant value for now, profile TBD)


    for (int i = 0; i < nlayer; i++) {// initialize alpha vals for all layers (constant value for now, profile TBD)
        alphas[i] = alpha;
    }
    for (int i = 0; i < nlayer; i++) {// initialize transmissivity vals for all layers (constant value for now, profile TBD)
        transmiss[i] = transmissivity;
    }

      /////////////////////////////////
     // Atmosphere initializations ///
    /////////////////////////////////


    // Initialize pressure levels
    for (int i = 0; i < nlayer + 1; i++){
        p_lev[i] = i * dp; // list of level pressures
    }

    //initialzie layer pressures
    for (int i = 0; i < nlayer; i++) {
        p_lay[i] = (p_lev[i] + p_lev[i+1]) / 2.0; // pressure for each layer as mean of the two pressure levels
    }

    // initialize temperature profile and geometric heights
    T_lev[nlayer] = T_lowest_init; // set initial temperature on the lowest level
    init_heights_temps(p_lev, T_lev, z_lev, GAMMA_DRY, nlayer+1); // z used as helper for temp init for now

   for (int i = 0; i < nlayer; i++){
    if (i == 0) {
        T_lay[i] = TtoTheta(p_lay[i], T_lowest_init, 1);
    }
    else {
        T_lay[i] = (T_lev[i] + T_lev[i+1]) / 2.0;
    }

    z_lay[i] = (z_lev[i] + z_lev[i+1]) / 2.0;
    }


    // calculate potential temperature
    for (int i = 0; i < nlayer; i++)
    {   
        theta_lay[i] = TtoTheta(p_lay[i], T_lay[i], 0); // calculate potential temperature for each layer given the starting temp an layer press
    }
 

    // Print results of initialization
    for (int i = 0; i < nlayer; i++)
    {
        printf("Layer %d: p = %.1f Pa, theta = %.2f K, temperature = %.2f K, altitude = %.2f m\n",
               i, p_lay[i], theta_lay[i], T_lay[i], z_lay[i]);
    }
 
      //////////////////////// 
     //// Main model loop ///
    ////////////////////////
    int count = 0;
    while (count < n_timesteps) {

        // do radiative transfer using only irradiance
        //irradiance_transfer(T_lay, alphas, transmiss, E_up, E_down, Eground, nlayer);

        // do radiative transfer using radiance
        radiance_transfer(T_lay, tau, mu, nmu, dmu, L_up, L_down, Eground,  nlayer);

        // adjust temperatures from level Energy transfer
        adjust_temps(L_down, L_up, L_net, L_abs, T_lay, theta_lay, p_lay, dt, dp, nlayer);
        
        // do convection
        convection(T_lay, theta_lay, p_lay, nlayer);


        // create CSV format output via print
        for (int i = 0; i < nlayer; i++)
        {

            printf("%i ; %.1f ; %.2f ; %.2f ; %.2f ; ",
                i, p_lay[i], theta_lay[i], T_lay[i], z_lay[i]);
            
            if (i == nlayer-1) {
                printf("\n");
            }
        }
        count++; // increment counter after each timestep
    }


    // Print results of model run
    printf("Model run completed after %i iterations \n", count);
    for (int i = 0; i < nlayer; i++)
    {
        printf("Layer %d: p = %.1f Pa, theta = %.2f K, temperature = %.2f K\n",
               i, p_lay[i], theta_lay[i], T_lay[i]);
    }


    return 0;
}