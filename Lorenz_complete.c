#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ===============================
// Parametri generali
// ===============================
#define NUM_STEPS 10000
#define T_STEPS 300
#define DELTAMAX 1.0
#define DELTA0 0.1
#define DT 0.01
#define TOL 0.4
#define SIGMA 10.0
#define BETA (8.0/3.0)

// Parametri rho
#define NUM_RHO 200
double RHO_START = 25.0;
double RHO_END   = 40.0;

// Parametri rumore
#define ADD_NOISE 1           // 1 = attivo, 0 = disattivo
#define NOISE_INTENSITY 0.01

// Parte deterministica azzerata
#define ZERO_DETERMINISTIC 0  // 1 = azzerata, 0 = attiva

// ===============================
// Funzioni helper
// ===============================
double rand_normal() {
    // Box-Muller per generare N(0,1)
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2*M_PI*u2);
}

void lorenz_system(double x[3], double rho, double dx[3]) {
    if (ZERO_DETERMINISTIC) {
        dx[0] = dx[1] = dx[2] = 0.0;
    } else {
        dx[0] = SIGMA * (x[1] - x[0]);
        dx[1] = x[0] * (rho - x[2]) - x[1];
        dx[2] = x[0] * x[1] - BETA * x[2];
    }
}

void wiener_increment(double dW[3]) {
    for(int i=0;i<3;i++) {
        dW[i] = ADD_NOISE ? NOISE_INTENSITY * sqrt(DT) * rand_normal() : 0.0;
    }
}

void lorenz_step(double x[3], double rho) {
    double dx[3], dW[3];
    lorenz_system(x, rho, dx);
    wiener_increment(dW);
    for(int i=0;i<3;i++)
        x[i] += dx[i]*DT + dW[i];
}

double jacobian_max(double x[3], double rho) {
    // Calcolo Lyapunov tradizionale semplificato
    double J[3][3];
    if (ZERO_DETERMINISTIC) {
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) J[i][j]=0.0;
    } else {
        J[0][0] = -SIGMA; J[0][1]=SIGMA; J[0][2]=0;
        J[1][0]=rho-x[2]; J[1][1]=-1;    J[1][2]=-x[0];
        J[2][0]=x[1];     J[2][1]=x[0];  J[2][2]=-BETA;
    }
    // Per semplicità restituiamo la norma massima della matrice
    double maxv = 0.0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            if(fabs(J[i][j])>maxv) maxv=fabs(J[i][j]);
    return maxv;
}

// ===============================
// Funzione principale
// ===============================
int main() {
    srand(time(NULL));

    // ===============================
    // Prepara array rho_values
    // ===============================
    double rho_values[NUM_RHO];
    for(int i=0;i<NUM_RHO;i++)
        rho_values[i] = RHO_START + i*(RHO_END - RHO_START)/(NUM_RHO-1);

    // ===============================
    // Loop su rho
    // ===============================
    for(int irho=0;irho<NUM_RHO;irho++){
        double rho = rho_values[irho];
        printf("--- rho = %.6f ---\n", rho);

        // Inizializza sistema
        double x[3] = {1.0,1.0,1.0};
        double v[3] = {1.0,0.0,0.0};
        double lyap_sum = 0.0;

        // ===============================
        // Genera traiettoria e calcola Lyapunov tradizionale
        // ===============================
        for(int step=0; step<NUM_STEPS; step++){
            lorenz_step(x, rho);

            // Lyapunov tradizionale
            double Jmax = jacobian_max(x, rho);
            for(int i=0;i<3;i++) v[i] = v[i] + DT*Jmax*v[i];
            double norm_v = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            for(int i=0;i<3;i++) v[i] /= norm_v;
            lyap_sum += log(norm_v);
        }

        double lyap_trad = lyap_sum / (NUM_STEPS*DT);
        double tempo_pred = INFINITY;

        if(ZERO_DETERMINISTIC){
            printf("Lyapunov tradizionale = %.6f | Tempo predicibilità = inf\n\n", lyap_trad);
        } else if(lyap_trad>0){
            tempo_pred = (1.0/lyap_trad)*log(DELTAMAX/DELTA0);
            printf("Lyapunov tradizionale = %.6f | Tempo predicibilità = %.3f\n\n", lyap_trad, tempo_pred);
        } else {
            printf("Lyapunov tradizionale = %.6f | Tempo predicibilità = inf\n\n", lyap_trad);
        }
    }

    return 0;
}
