#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "funcs.h"

//menu 1
static double g_R  = 1000.0;     
static double g_C  = 1e-6;       
static double g_Vin = 5.0;       
static double g_t_max = 0.01;   
static double g_dt    = 0.0001;

//menu 2
extern double g_R;
extern double g_C;
extern double g_Vin;
extern double g_t_max;
extern double g_dt;

 //menu_item_1: Set RC circuit parameters
void menu_item_1(void) {
   double R, C, Vin, tmax, dt;   // Unit 1: local variables (automatic storage)

    printf("\n=== Set Circuit Parameters ===\n");

    // Unit 1: scanf + & (address operator) 
    printf("Enter resistance R (ohms): ");
    if (scanf("%lf", &R) != 1) return;

    printf("Enter capacitance C (farads): ");
    if (scanf("%lf", &C) != 1) return;

    printf("Enter input voltage Vin (volts): ");
    if (scanf("%lf", &Vin) != 1) return;

    printf("Enter simulation max time t_max (seconds): ");
    if (scanf("%lf", &tmax) != 1) return;

    printf("Enter time step dt (seconds): ");
    if (scanf("%lf", &dt) != 1) return;

    // Unit 1: conditional checks 
    if (R <= 0 || C <= 0 || tmax <= 0 || dt <= 0) {
        printf("\nError: All parameters must be positive numbers.\n");
        return;
    }


    g_R = R;
    g_C = C;
    g_Vin = Vin;
    g_t_max = tmax;
    g_dt = dt;

    printf("\n>>> Updated parameters:\n");
    printf("R      = %g ohms\n", g_R);
    printf("C      = %g F\n",   g_C);
    printf("Vin    = %g V\n",   g_Vin);
    printf("t_max  = %g s\n",   g_t_max);
    printf("dt     = %g s\n",   g_dt);

    printf("\nParameters successfully updated.\n");
}

//menu_item_2: Simulate RC Step Response (0 → Vin)
void menu_item_2(void) {
    printf("\n=== Step Response: 0 → Vin ===\n");

    /* Unit 1: use local variables */
    double R = g_R;
    double C = g_C;
    double Vin = g_Vin;
    double t_max = g_t_max;
    double dt = g_dt;

    double tau = R * C;   // Time constant τ

    printf("Using R = %g ohms, C = %g F, Vin = %g V\n", R, C, Vin);
    printf("Time constant τ = R*C = %g s\n", tau);
    printf("Simulating 0 → %.6f seconds with dt = %.6f seconds...\n",
           t_max, dt);

    
    int steps = (int)(t_max / dt) + 1;

    double voltage[steps];   // Variable Length Array (C99)

    // Unit 1: for loop 
    for (int i = 0; i < steps; i++) {
        double t = i * dt;

        // Unit 1: math.h exp() 
        voltage[i] = Vin * (1.0 - exp(-t / tau));
    }

    printf("\n t (s)\t\tVc(t) (V)\n");
    printf("-----------------------------\n");

    /* Print every 5 points */
    for (int i = 0; i < steps; i += 5) {
        double t = i * dt;
        printf("%10.6f\t%10.6f\n", t, voltage[i]);
    }

    printf("\nSimulation complete.\n");
}

//menu_item_3: Discharge response (V0 → 0)
void menu_item_3(void) {
    printf("\n=== Discharge Response: V0 → 0 ===\n");

    
    double R     = g_R;
    double C     = g_C;
    double t_max = g_t_max;
    double dt    = g_dt;

    if (R <= 0 || C <= 0) {
        printf("Error: R and C must be set to positive values first (use menu 1).\n");
        return;
    }

    double V0;
    printf("Enter initial capacitor voltage V0 (volts): ");
    if (scanf("%lf", &V0) != 1) {
        printf("Input error.\n");
        return;
    }

    if (V0 < 0) {
        printf("Warning: V0 is negative, using its absolute value.\n");
        V0 = -V0;
    }

    double tau = R * C;   /* time constant */

    printf("\nUsing R = %g ohms, C = %g F\n", R, C);
    printf("Time constant τ = R*C = %g s\n", tau);
    printf("Simulating discharge from V0 = %g V over 0 → %.6f s with dt = %.6f s\n",
           V0, t_max, dt);

    /* Unit 1: array + pointer
       Compute number of steps
    */
    int steps = (int)(t_max / dt) + 1;
    if (steps <= 0) {
        printf("Error: invalid simulation settings (t_max or dt).\n");
        return;
    }

    double voltage[steps];   // Array as a pointer to double 

    /* Unit 1: for loop + exp() */
    for (int i = 0; i < steps; i++) {
        double t = i * dt;
        voltage[i] = V0 * exp(-t / tau);
    }

    printf("\n t (s)\t\tVc(t) (V)\n");
    printf("-----------------------------\n");

    // Unit 1: pointer
    double *p = voltage;

    //Print every 5 points 
    for (int i = 0; i < steps; i += 5) {
        double t = i * dt;
        double vc = *(p + i);  
        printf("%10.6f\t%10.6f\n", t, vc);
    }

    printf("\nDischarge simulation complete.\n");
}

//menu_item_4: Export step/discharge data to CSV
void menu_item_4(void) {
    
    char filename[100];

    printf("\n=== Export BOTH Step & Discharge Responses to CSV ===\n");
    printf("Enter output filename (e.g. both.csv): ");

    if (scanf("%99s", filename) != 1) {
        printf("Input error.\n");
        return;
    }

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: could not open file '%s' for writing.\n", filename);
        return;
    }

    double R     = g_R;
    double C     = g_C;
    double Vin   = g_Vin;
    double t_max = g_t_max;
    double dt    = g_dt;

    if (R <= 0 || C <= 0 || t_max <= 0 || dt <= 0) {
        printf("Error: invalid circuit or simulation parameters.\n");
        fclose(fp);
        return;
    }

    double tau = R * C;

    /* Assume discharge starts from V0 = Vin */
    double V0 = Vin;

    int steps = (int)(t_max / dt) + 1;
    if (steps <= 0) {
        printf("Error: invalid number of steps.\n");
        fclose(fp);
        return;
    }

    /* CSV header */
    fprintf(fp, "t,step_vc,discharge_vc\n");

    /* Compute both responses */
    for (int i = 0; i < steps; i++) {
        double t = i * dt;

        double vc_step = Vin * (1.0 - exp(-t / tau));   // Step response
        double vc_dis  = V0  * exp(-t / tau);           // Discharge response

        fprintf(fp, "%.6f,%.6f,%.6f\n", t, vc_step, vc_dis);
    }

    fclose(fp);

    printf("Export completed.\n");
    printf("File '%s' now contains BOTH curves (step & discharge).\n", filename);
    printf("You can open it in Excel or MATLAB to plot both curves.\n");
}

//menu_item_5: Analyse time constant τ
void menu_item_5(void) {
     printf("\n=== Analyse Time Constant τ ===\n");

    double R     = g_R;
    double C     = g_C;
    double Vin   = g_Vin;
    double t_max = g_t_max;
    double dt    = g_dt;

    /* Basic parameter checks */
    if (R <= 0 || C <= 0) {
        printf("Error: R and C must be positive. Please set them in menu 1.\n");
        return;
    }
    if (t_max <= 0 || dt <= 0) {
        printf("Error: simulation settings (t_max, dt) must be positive.\n");
        return;
    }
    if (Vin <= 0) {
        printf("Warning: Vin <= 0, step response will not rise. Set Vin > 0 in menu 1.\n");
        return;
    }

    /* Theoretical time constant */
    double tau_theory = R * C;

    printf("Current circuit parameters:\n");
    printf("R      = %g ohms\n", R);
    printf("C      = %g F\n",   C);
    printf("Vin    = %g V\n",   Vin);
    printf("t_max  = %g s, dt = %g s\n", t_max, dt);
    printf("\nTheoretical time constant τ_theory = R * C = %g s\n", tau_theory);

    /* Target voltage: 0.632 Vin */
    const double target_ratio = 1.0 - exp(-1.0);
    double target_v = target_ratio * Vin;

    int steps = (int)(t_max / dt) + 1;
    if (steps <= 1) {
        printf("Error: simulation time window is too short.\n");
        return;
    }

    /* Numerical step response: find where Vc ≈ 0.632 Vin */
    double prev_t  = 0.0;
    double prev_vc = 0.0;
    int found = 0;
    double t_est = 0.0;

    for (int i = 1; i < steps; i++) {
        double t  = i * dt;
        double vc = Vin * (1.0 - exp(-t / tau_theory));

        /* Linear interpolation between two time points */
        if (vc >= target_v) {
            found = 1;
            if (vc == prev_vc) {
                t_est = t;
            } else {
                t_est = prev_t + (target_v - prev_vc) * (t - prev_t) / (vc - prev_vc);
            }
            break;
        }

        prev_t  = t;
        prev_vc = vc;
    }

    if (!found) {
        printf("\nWarning: Within t_max = %g s, Vc(t) did not reach 0.632*Vin.\n", t_max);
        printf("Increase t_max (e.g. ~3*τ_theory) in menu 1 or 7 and try again.\n");
        return;
    }

    double tau_est = t_est;

    printf("\nEstimated time constant from step response:\n");
    printf("τ_est ≈ %g s (time at which Vc ≈ 0.632 * Vin)\n", tau_est);

    double error = 0.0;
    if (tau_theory != 0.0) {
        error = (tau_est - tau_theory) / tau_theory * 100.0;
    }

    printf("\nComparison:\n");
    printf("τ_theory = %g s\n", tau_theory);
    printf("τ_est    = %g s\n", tau_est);
    printf("Relative error ≈ %.4f %%\n", error);
}

//menu_item_6: Frequency response at a given frequency
void menu_item_6(void) {
     printf("\n=== Frequency Response of RC Circuit ===\n");

    double R = g_R;
    double C = g_C;

    if (R <= 0 || C <= 0) {
        printf("Error: R and C must be positive. Please set them in menu 1.\n");
        return;
    }

    double RC = R * C;
    printf("Current circuit parameters:\n");
    printf("R = %g ohms, C = %g F, RC = %g s\n", R, C, RC);
    printf("\nEnter frequency values in Hz (f > 0).\n");
    printf("Enter 0 or a negative value to return to the main menu.\n\n");

    while (1) {
        double f;
        printf("Frequency f (Hz): ");

        if (scanf("%lf", &f) != 1) {
            printf("Input error. Stopping frequency response analysis.\n");
            return;
        }

        if (f <= 0.0) {
            printf("Exiting frequency response analysis.\n");
            break;
        }

        double omega = 2.0 * M_PI * f;

        double omegaRC = omega * RC;
        double mag = 1.0 / sqrt(1.0 + omegaRC * omegaRC);

        double gain_dB = 20.0 * log10(mag);

        double phase_rad = -atan(omegaRC);
        double phase_deg = phase_rad * 180.0 / M_PI;

        printf("\nAt f = %g Hz:\n", f);
        printf("  |H(jω)|   = %g\n", mag);
        printf("  Gain(dB)  = %g dB\n", gain_dB);
        printf("  Phase     = %g degrees\n\n", phase_deg);
    }
}

//menu_item_7: Change simulation settings (t_max, dt)
void menu_item_7(void) {
    printf("\n=== Change Simulation Settings (t_max, dt) ===\n");

    printf("Current settings:\n");
    printf("  t_max = %g s\n", g_t_max);
    printf("  dt    = %g s\n", g_dt);

    double new_tmax, new_dt;

    printf("\nEnter new t_max (seconds): ");
    if (scanf("%lf", &new_tmax) != 1) {
        printf("Input error.\n");
        return;
    }

    printf("Enter new dt (seconds): ");
    if (scanf("%lf", &new_dt) != 1) {
        printf("Input error.\n");
        return;
    }

    /* Basic checks */
    if (new_tmax <= 0.0 || new_dt <= 0.0) {
        printf("Error: t_max and dt must both be positive.\n");
        return;
    }

    if (new_dt >= new_tmax) {
        printf("Error: dt should be smaller than t_max to have multiple steps.\n");
        return;
    }

    int steps = (int)(new_tmax / new_dt) + 1;

    printf("\nWith these settings you will have about %d simulation points.\n", steps);

    if (steps > 200000) {
        printf("Warning: This is a very large number of points.\n");
        printf("The simulation and CSV export may be slow or generate very large files.\n");
    } else if (steps < 10) {
        printf("Warning: Very few points — you may not see a smooth curve.\n");
    }

    /* Update global variables */
    g_t_max = new_tmax;
    g_dt    = new_dt;

    printf("\nSimulation settings updated:\n");
    printf("  t_max = %g s\n", g_t_max);
    printf("  dt    = %g s\n", g_dt);
}
