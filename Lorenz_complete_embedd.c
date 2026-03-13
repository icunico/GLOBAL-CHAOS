#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_STEPS 10000
#define T_STEPS 300
#define DELTAMAX 2.0
#define DELTA0 0.1
#define DT 0.01
#define TOL 0.4

#define SIGMA 10.0
#define BETA (8.0/3.0)

// ===== Parametri rho =====
#define RHO_MIN 0.0
#define RHO_MAX 40.0
#define RHO_STEP 0.5

// ===== Flags =====
#define ADD_NOISE 0
#define ZERO_DETERMINISTIC 0
#define USE_EMBEDDING 1  // 0 = originale, 1 = embedding 3D

double noise_intensity = 0.01;

// ===== Structs =====
typedef struct { double x[3]; } Vec3;

typedef struct {
    Vec3 *x_vals;
    Vec3 *y_vals;
    int length;
} RecurrencePair;

// ===== Funzioni utili =====
double randn() {
    double u1 = ((double)rand()) / RAND_MAX;
    double u2 = ((double)rand()) / RAND_MAX;
    return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

double norm(Vec3 a, Vec3 b) {
    double s = 0;
    for(int i=0;i<3;i++){
        double d = a.x[i] - b.x[i];
        s += d*d;
    }
    return sqrt(s);
}

// ===== Sistema di Lorenz =====
Vec3 lorenz_system(Vec3 x, double rho) {
    Vec3 out;
    if(ZERO_DETERMINISTIC){
        out.x[0]=0; out.x[1]=0; out.x[2]=0;
        return out;
    }
    out.x[0] = SIGMA*(x.x[1]-x.x[0]);
    out.x[1] = x.x[0]*(rho - x.x[2]) - x.x[1];
    out.x[2] = x.x[0]*x.x[1] - BETA*x.x[2];
    return out;
}

Vec3 wiener_increment() {
    Vec3 w;
    for(int i=0;i<3;i++)
        w.x[i] = noise_intensity*sqrt(DT)*randn();
    return w;
}

Vec3 lorenz_system_sde(Vec3 x, double rho){
    Vec3 drift = lorenz_system(x,rho);
    Vec3 diffusion;
    if(ADD_NOISE) diffusion = wiener_increment();
    else diffusion.x[0]=diffusion.x[1]=diffusion.x[2]=0;

    Vec3 out;
    for(int i=0;i<3;i++)
        out.x[i] = x.x[i] + drift.x[i]*DT + diffusion.x[i];
    return out;
}

// ===== Jacobiano =====
void jacobian(Vec3 x, double rho, double J[3][3]){
    if(ZERO_DETERMINISTIC){
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                J[i][j]=0;
        return;
    }
    J[0][0]=-SIGMA; J[0][1]=SIGMA; J[0][2]=0;
    J[1][0]=rho - x.x[2]; J[1][1]=-1; J[1][2]=-x.x[0];
    J[2][0]=x.x[1]; J[2][1]=x.x[0]; J[2][2]=-BETA;
}

// ===== Ricorrenze con embedding corretto =====
RecurrencePair* select_y_by_recurrence(Vec3 *series, int N, int *n_pairs){
    double soglia_min = DELTA0*(1-TOL);
    double soglia_max = DELTA0*(1+TOL);
    int max_rec = 500;
    RecurrencePair *pairs = malloc(max_rec*sizeof(RecurrencePair));
    int rec_found = 0;
    int last_t1_used = -T_STEPS;

    int emb_dim = 3;          // dimensione embedding
    int tau_emb = 10;          // ritardo embedding
    int max_index = N - (USE_EMBEDDING ? (emb_dim-1)*tau_emb : 0);

    for(int t1=0; t1<max_index-T_STEPS; t1++){
        if(t1-last_t1_used < T_STEPS) continue;
        for(int t2=t1+T_STEPS; t2<max_index; t2++){
            Vec3 v1 = {0}, v2 = {0};
            if(USE_EMBEDDING){
                // Creazione embedding solo lungo x
                for(int k=0;k<emb_dim;k++){
                    v1.x[k] = series[t1 + k*tau_emb].x[0];
                    v2.x[k] = series[t2 + k*tau_emb].x[0];
                }
            } else {
                v1 = series[t1];
                v2 = series[t2];
            }

            double d = norm(v1, v2);
            if(d>=soglia_min && d<=soglia_max){
                int max_len = max_index - (t1>t2 ? t1 : t2);
                pairs[rec_found].x_vals = malloc(max_len*sizeof(Vec3));
                pairs[rec_found].y_vals = malloc(max_len*sizeof(Vec3));
                for(int k=0;k<max_len;k++){
                    if(USE_EMBEDDING){
                        // Copia solo x, y,z=0
                        pairs[rec_found].x_vals[k].x[0] = series[t1+k].x[0];
                        pairs[rec_found].x_vals[k].x[1] = 0;
                        pairs[rec_found].x_vals[k].x[2] = 0;
                        pairs[rec_found].y_vals[k].x[0] = series[t2+k].x[0];
                        pairs[rec_found].y_vals[k].x[1] = 0;
                        pairs[rec_found].y_vals[k].x[2] = 0;
                    } else {
                        pairs[rec_found].x_vals[k] = series[t1+k];
                        pairs[rec_found].y_vals[k] = series[t2+k];
                    }
                }
                pairs[rec_found].length = max_len;
                rec_found++;
                last_t1_used = t1;
                if(rec_found>=max_rec) break;
                break;
            }
        }
        if(rec_found>=max_rec) break;
    }
    *n_pairs = rec_found;
    return pairs;
}

// ===== Calcolo Lyapunov tau =====
void lyap_tau_continuous_multi(RecurrencePair *pairs,int n_pairs,double *Tpred,double *Tcrit){
    int tau_count=0;
    double tau_sum=0,d0_sum=0;
    double soglia_min = DELTA0*(1-TOL);
    double soglia_max = DELTA0*(1+TOL);

    for(int p=0;p<n_pairs;p++){
        int counting=0,tau_current=0;
        double d0_current=0;
        for(int i=0;i<pairs[p].length;i++){
            double d = norm(pairs[p].x_vals[i],pairs[p].y_vals[i]);
            if(!counting){
                if(d>=soglia_min && d<=soglia_max){
                    counting=1; tau_current=0; d0_current=d;
                }
            } else {
                tau_current++;
                if(d>=DELTAMAX){
                    tau_sum += tau_current;
                    d0_sum += d0_current;
                    tau_count++;
                    counting=0;
                }
            }
        }
    }

    if(tau_count==0){ *Tpred = NAN; *Tcrit = NAN; return; }

    double mean_tau = tau_sum/tau_count;
    double mean_d0  = d0_sum/tau_count;
    *Tcrit = mean_tau*DT;
    *Tpred = (*Tcrit) * log(DELTAMAX/DELTA0)/log(DELTAMAX/mean_d0);
}

// ===== MAIN =====
int main(){
    srand(time(NULL));

    int rho_count = (int)((RHO_MAX-RHO_MIN)/RHO_STEP + 1);
    double *rho_values = malloc(rho_count*sizeof(double));
    for(int i=0;i<rho_count;i++)
        rho_values[i] = RHO_MIN + i*RHO_STEP;

    double *tempo_pred_trad_list = malloc(rho_count*sizeof(double));
    double *tempo_pred_tau_list  = malloc(rho_count*sizeof(double));
    double *tempo_critico_list   = malloc(rho_count*sizeof(double));

    for(int i=0;i<rho_count;i++){
        double rho = rho_values[i];
        printf("\n--- Simulazione %d/%d: rho=%.2f ---\n",i+1,rho_count,rho);

        Vec3 *series = malloc((NUM_STEPS+1)*sizeof(Vec3));
        series[0].x[0]=1; series[0].x[1]=1; series[0].x[2]=1;

        for(int step=0;step<NUM_STEPS;step++)
            series[step+1] = lorenz_system_sde(series[step],rho);

        // Lyapunov tradizionale
        double v[3]={1,0,0};
        double lyap_sum=0;
        for(int step=0;step<NUM_STEPS;step++){
            double J[3][3];
            jacobian(series[step],rho,J);
            double newv[3];
            for(int k=0;k<3;k++)
                newv[k] = v[k] + DT*(v[0]*J[k][0] + v[1]*J[k][1] + v[2]*J[k][2]);
            double normv = sqrt(newv[0]*newv[0]+newv[1]*newv[1]+newv[2]*newv[2]);
            for(int k=0;k<3;k++) v[k]=newv[k]/normv;
            lyap_sum += log(normv);
        }
        double lyap_trad = lyap_sum/(NUM_STEPS*DT);
        double tempo_pred_trad = (ZERO_DETERMINISTIC || lyap_trad<=0) ? INFINITY : (1.0/lyap_trad)*log(DELTAMAX/DELTA0);
        tempo_pred_trad_list[i] = tempo_pred_trad;

        // Ricorrenze e tau
        int n_pairs;
        RecurrencePair *pairs = select_y_by_recurrence(series,NUM_STEPS,&n_pairs);
        double Tpred_tau,Tcrit;
        lyap_tau_continuous_multi(pairs,n_pairs,&Tpred_tau,&Tcrit);
        tempo_pred_tau_list[i] = Tpred_tau;
        tempo_critico_list[i]  = Tcrit;

        printf("Lyapunov tradizionale: %f\n",lyap_trad);
        printf("Tempo predicibilita tradizionale: %f\n",tempo_pred_trad);
        printf("Tempo predicibilita tau: %f\n",Tpred_tau);
        printf("Tempo critico: %f\n",Tcrit);

        for(int p=0;p<n_pairs;p++){
            free(pairs[p].x_vals);
            free(pairs[p].y_vals);
        }
        free(pairs);
        free(series);
    }

    FILE *f = fopen("tempo_predicibilita_vs_rho.dat","w");
    fprintf(f,"# rho tempo_pred_trad tempo_pred_tau tempo_critico\n");
    for(int i=0;i<rho_count;i++)
        fprintf(f,"%f %f %f %f\n", rho_values[i],tempo_pred_trad_list[i],tempo_pred_tau_list[i],tempo_critico_list[i]);
    fclose(f);

    free(tempo_pred_trad_list);
    free(tempo_pred_tau_list);
    free(tempo_critico_list);
    free(rho_values);

    printf("\n=== SIMULAZIONE COMPLETATA ===\n");
    return 0;
}
