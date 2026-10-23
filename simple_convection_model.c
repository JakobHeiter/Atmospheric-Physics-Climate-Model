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

  ////////////////////////////
 /// Main program function //
////////////////////////////


int main(int argc, char **argv)
{
    // Model constraints
    int nlayer = 10;
    int n_timesteps = 10000;

    // Define atmosphere TODO: sinnvoll sortieren!
    double Eground = 240; // default for incoming energy on ground (W/m^2)
    double dt = 3600; // time increment for model (s) 
    double dp = 100000.0 / nlayer; // pressure difference betw. levels (Pa)
    double T_lowest_init = 288; // initial temperature of lowest layer (K)

    double alpha = 0.95; // absorption coefficient (const for all layers currently)
    double transmissivity = 1-alpha;

    double p_lev[nlayer + 1]; // pressure levels (Pa)
    double T_lev[nlayer + 1];   // level temperatures (K)
    double z_lev[nlayer + 1]; // geometric height levels (m)

    double theta_lay[nlayer]; // layer potential temperatures
    double T_lay[nlayer]; // layer temperatures (K)
    double p_lay[nlayer]; // mean assumed pressure for each layer (Pa)
    double z_lay[nlayer]; // geometric height of each layer (m)


    double E_up [nlayer+1]; // upward energy flux at each level (W/m^2)
    double E_down [nlayer+1]; // downward energy flux at each level (W/m^2)
    double E_net [nlayer+1]; // net energy flux at each level(W/m^2)

    double E_abs [nlayer]; // LAYER QTY!! Net absorbed energy at each layer (E_net[i] - E_net[i+1])

    double alphas [nlayer]; // absorption coefficients of each layer
    double transmiss [nlayer]; // transmissivity of each layer (t = 1-alpha)

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
    for (int i = 0; i < nlayer + 1; i++) //
    {
        p_lev[i] = i * dp; // list of level pressures
    }

    //initialzie layer pressures
    for (int i = 0; i < nlayer; i++) 
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

        // Calculate the heating rate and heat lowest layer
        //double dT = heating_rate(Eground, dt, dp);
        //T_lay[nlayer - 1] += dT;// Update the temperature of the lowermost layer
        //theta_lay[nlayer - 1] = TtoTheta(p_lay[nlayer - 1], T_lay[nlayer - 1], 0);// Update the potential temperature of the lowermost layer


        // Calculate Flux boundary conditions
        E_up[nlayer] = Eground; // flux from lowest layer is just ground flux (+ emission of layer itself which is calculated in the upward flux loop)//SIGMA * pow(T_lay[nlayer-1], 4); // lowest layer heating the layer above (Boltzmann)
        E_down[0] = 0; // highest layer does not receive significant amounts of energy from space

        // calculate upwards fluxes at each level
        for (int i = nlayer-1; i > -1 ; i--){ // TODO: check correct loop boundaries and indicecs
            E_up[i] = E_up[i+1] * transmiss[i] + SIGMA * pow(T_lay[i], 4) * alphas[i];
        }

        //calculate downwards fluxes at each level
        for (int i = 0; i < nlayer; i++){ // so, anstatt i < nlayer +1, weil wir vom oberstem level bis zum boden E_down berechnen
            E_down[i+1] =  E_down[i] * transmiss[i] + SIGMA* pow(T_lay[i], 4) * alphas[i];
        }

        //calculate E_net for each level (E_net = E_dn - E_up)
        for (int i = 0; i < nlayer + 1;  i++){
            E_net[i] = E_down[i] - E_up[i];
        }

        //calculate E_abs for each layer (E_abs = Enet von oben - E_net von unten) problematischster Teil
        for (int i = 0; i < nlayer; i++){   // Layer quantity
            E_abs[i] = E_net[i] - E_net[i+1]; // minus makes sense as fluxes pointing upwards are negative --> - (a negative quantity) > 0, basically adding energy
            //if (i == nlayer - 1){
            //    E_abs[i] = Eground - E_net[i + 1]; // TODO: ground heating included here! Must remove above!
            //} this is all removed bc at lowest layer E_up[nlayer] = Eground; i.e. surface flux is already accounted for!!
        }

        // update temperature for each layer using E_abs 
        for (int i = 0; i < nlayer; i++){
            double dT = heating_rate(E_abs[i], dt, dp);
            T_lay[i] += dT;
        }
        
        // update potential temperatures for each layer using E_abs
        for (int i = 0; i < nlayer; i++){
            theta_lay[i] = TtoTheta(p_lay[i], T_lay[i], 0);// Update the potential temperature of the layers
        }
        
        // do convection
        convection(T_lay, theta_lay, p_lay, nlayer);
        

        // create CSV format output via print
        for (int i = 0; i < nlayer; i++)
        {

            printf("%i ; %.1f ; %.2f ; %.2f ;",
                i, p_lay[i], theta_lay[i], T_lay[i]);
            
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