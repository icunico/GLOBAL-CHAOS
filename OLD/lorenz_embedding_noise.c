#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// ===============================
// Parametri globali
// ===============================
#define NUM_STEPS 100000
#define T_STEPS 300
#define DELTA_MAX 2.0
#define DELTA0 0.1
#define DT 0.01
#define TOL 0.4
#define SIGMA 10.0
#define BETA (8.0/3.0)
#define NUM_RHO 10
#define MAX_RECURRENCES 500

// Parametri rumore
#define ADD_NOISE 1  // 1 per aggiungere rumore, 0 per non aggiungere
#define NOISE_INTENSITY 1.
#define USE_EMBEDDING 0 // 1 per embedding, 0 per output diretto

#ifndef NAN
#define NAN (0.0/0.0)
#endif

// ===============================
// Strutture dati
// ===============================
typedef struct {
    double x, y, z;
} Point3D;

typedef struct {
    Point3D *points;
    int length;
} Trajectory;

typedef struct {
    double **matrix;
    int rows, cols;
} Matrix;

typedef struct {
    Trajectory x_traj;
    Trajectory y_traj;
} RecurrencePair;

typedef struct {
    RecurrencePair *pairs;
    int count;
} RecurrenceList;

// ===============================
// Generatore di numeri casuali (Box-Muller)
// ===============================
double randn() {
    static int has_spare = 0;
    static double spare;
    
    if (has_spare) {
        has_spare = 0;
        return spare;
    }
    
    has_spare = 1;
    static double u, v, mag;
    do {
        u = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        v = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        mag = u * u + v * v;
    } while (mag >= 1.0 || mag == 0.0);
    
    mag = sqrt(-2.0 * log(mag) / mag);
    spare = v * mag;
    return u * mag;
}

// ===============================
// Funzioni ausiliarie
// ===============================
double norm3d(Point3D p1, Point3D p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

Point3D add_points(Point3D p1, Point3D p2) {
    Point3D result = {p1.x + p2.x, p1.y + p2.y, p1.z + p2.z};
    return result;
}

Point3D scale_point(Point3D p, double scale) {
    Point3D result = {p.x * scale, p.y * scale, p.z * scale};
    return result;
}

double vector_norm(Point3D p) {
    return sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
}

Point3D normalize_vector(Point3D p) {
    double norm = vector_norm(p);
    if (norm > 1e-10) {
        return scale_point(p, 1.0/norm);
    }
    return p;
}

Matrix* create_matrix(int rows, int cols) {
    Matrix *m = malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    m->matrix = malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        m->matrix[i] = malloc(cols * sizeof(double));
    }
    return m;
}

void free_matrix(Matrix *m) {
    if (m) {
        for (int i = 0; i < m->rows; i++) {
            free(m->matrix[i]);
        }
        free(m->matrix);
        free(m);
    }
}

Point3D matrix_vector_mult(Matrix *m, Point3D v) {
    Point3D result;
    result.x = m->matrix[0][0]*v.x + m->matrix[0][1]*v.y + m->matrix[0][2]*v.z;
    result.y = m->matrix[1][0]*v.x + m->matrix[1][1]*v.y + m->matrix[1][2]*v.z;
    result.z = m->matrix[2][0]*v.x + m->matrix[2][1]*v.y + m->matrix[2][2]*v.z;
    return result;
}

void free_trajectory(Trajectory *traj) {
    if (traj && traj->points) {
        free(traj->points);
        traj->points = NULL;
        traj->length = 0;
    }
}

void free_recurrence_list(RecurrenceList *rec_list) {
    if (rec_list && rec_list->pairs) {
        for (int i = 0; i < rec_list->count; i++) {
            free_trajectory(&rec_list->pairs[i].x_traj);
            free_trajectory(&rec_list->pairs[i].y_traj);
        }
        free(rec_list->pairs);
        rec_list->pairs = NULL;
        rec_list->count = 0;
    }
}

// ===============================
// Sistema di Lorenz con rumore
// ===============================
Point3D lorenz_system(Point3D x, double rho) {
    Point3D result;
    result.x = 0.0*(SIGMA * (x.y - x.x));  // dX/dt
    result.y = 0.0*(x.x * (rho - x.z) - x.y);  // dY/dt
    result.z = 0.0*(x.x * x.y - BETA * x.z);  // dZ/dt
    return result;
}

Point3D wiener_increment(double dt, double noise_intensity) {
    Point3D result;
    double sqrt_dt = sqrt(dt);
    result.x = noise_intensity * sqrt_dt * randn();
    result.y = noise_intensity * sqrt_dt * randn();
    result.z = noise_intensity * sqrt_dt * randn();
    return result;
}

Point3D lorenz_system_sde(Point3D x, double rho, double dt, double noise_intensity) {
    // Termine deterministico (drift)
    Point3D drift = lorenz_system(x, rho);
    
    // Termine stocastico (diffusione) - processo di Wiener
    Point3D diffusion;
    if (ADD_NOISE) {
        diffusion = wiener_increment(dt, noise_intensity);
    } else {
        diffusion.x = diffusion.y = diffusion.z = 0.0;
    }
    
    // Integrazione secondo Itō: dx = f(x)dt + g(x)dW
    Point3D x_new = add_points(x, add_points(scale_point(drift, dt), diffusion));
    return x_new;
}

Trajectory create_3d_embedding(double *trajectory_1d, int length, int tau, int embedding_dim) {
    Trajectory embedded;
    embedded.length = length - (embedding_dim - 1) * tau;
    
    if (embedded.length <= 0) {
        printf("Errore: Tau troppo grande per l'embedding\n");
        embedded.points = NULL;
        embedded.length = 0;
        return embedded;
    }
    
    embedded.points = malloc(embedded.length * sizeof(Point3D));
    
    for (int i = 0; i < embedded.length; i++) {
        embedded.points[i].x = trajectory_1d[i];
        embedded.points[i].y = trajectory_1d[i + tau];
        embedded.points[i].z = trajectory_1d[i + 2 * tau];
    }
    
    printf("  Embedding creato: %d -> %d punti\n", length, embedded.length);
    printf("  Parametri: tau=%d, dim=%d\n", tau, embedding_dim);
    
    return embedded;
}

Matrix* jacobian(Point3D x, double rho) {
    Matrix *J = create_matrix(3, 3);
    
    J->matrix[0][0] = -SIGMA;
    J->matrix[0][1] = SIGMA;
    J->matrix[0][2] = 0;
    
    J->matrix[1][0] = rho - x.z;
    J->matrix[1][1] = -1;
    J->matrix[1][2] = -x.x;
    
    J->matrix[2][0] = x.y;
    J->matrix[2][1] = x.x;
    J->matrix[2][2] = -BETA;
    
    return J;
}

// ===============================
// Selezione ricorrenze
// ===============================
RecurrenceList select_y_by_recurrence(Trajectory *traj, double delta0, double tol, 
                                     int t_steps, int max_recurrences) {
    printf("  Cercando ricorrenze nella traiettoria di %d punti...\n", traj->length);
    
    RecurrenceList rec_list = {NULL, 0};
    rec_list.pairs = malloc(max_recurrences * sizeof(RecurrencePair));
    
    int recurrences_found = 0;
    double soglia_min = delta0 * (1 - tol);
    double soglia_max = delta0 * (1 + tol);
    printf("  Soglia ricorrenze: %.4f <= dist <= %.4f\n", soglia_min, soglia_max);
    
    int last_t1_used = -t_steps;
    
    for (int t1 = 0; t1 < traj->length - t_steps && recurrences_found < max_recurrences; t1++) {
        // Verifica che t1 sia almeno t_steps passi dopo l'ultima ricorrenza trovata
        if (t1 - last_t1_used < t_steps) {
            continue;
        }
        
        for (int t2 = t1 + t_steps; t2 < traj->length; t2++) {
            double dist = norm3d(traj->points[t2], traj->points[t1]);
            
            if (dist >= soglia_min && dist <= soglia_max) {
                int max_len = traj->length - (t1 > t2 ? t1 : t2);
                last_t1_used = t1;
                
                // Crea nuove traiettorie per questa ricorrenza
                rec_list.pairs[recurrences_found].x_traj.length = max_len;
                rec_list.pairs[recurrences_found].y_traj.length = max_len;
                rec_list.pairs[recurrences_found].x_traj.points = malloc(max_len * sizeof(Point3D));
                rec_list.pairs[recurrences_found].y_traj.points = malloc(max_len * sizeof(Point3D));
                
                for (int i = 0; i < max_len; i++) {
                    rec_list.pairs[recurrences_found].x_traj.points[i] = traj->points[t1 + i];
                    rec_list.pairs[recurrences_found].y_traj.points[i] = traj->points[t2 + i];
                }
                
                recurrences_found++;
                printf("  Ricorrenza %d trovata: t1=%d, t2=%d, dist=%.4f, lunghezza=%d\n", 
                       recurrences_found, t1, t2, dist, max_len);
                
                if (recurrences_found >= max_recurrences) {
                    printf("  Raggiunto limite massimo di ricorrenze (%d)\n", max_recurrences);
                    break;
                }
                break;
            }
        }
    }
    
    rec_list.count = recurrences_found;
    
    if (recurrences_found == 0) {
        printf("  Nessuna ricorrenza trovata!\n");
        free(rec_list.pairs);
        rec_list.pairs = NULL;
    } else {
        printf("  Totale ricorrenze trovate: %d\n", rec_list.count);
    }
    
    return rec_list;
}

// ===============================
// Lyapunov tau multi-ricorrenze
// ===============================
double lyap_tau_continuous_multi(RecurrenceList *rec_list, double delta_max, 
                                double delta0, double tol) {
    if (!rec_list || rec_list->count == 0) {
        printf("  Nessuna ricorrenza fornita per il calcolo!\n");
        return NAN;
    }
    
    printf("  Calcolo Lyapunov tau su %d ricorrenze...\n", rec_list->count);
    
    double *all_tau_list = malloc(10000 * sizeof(double));
    double *all_d0_list = malloc(10000 * sizeof(double));
    int total_tau_count = 0;
    
    double soglia_min = delta0 * (1 - tol);
    double soglia_max = delta0 * (1 + tol);
    
    // Processa ogni coppia di ricorrenze
    for (int pair_idx = 0; pair_idx < rec_list->count; pair_idx++) {
        Trajectory *x_vals = &rec_list->pairs[pair_idx].x_traj;
        Trajectory *y_vals = &rec_list->pairs[pair_idx].y_traj;
        
        printf("    Processando ricorrenza %d/%d (%d punti)...\n", 
               pair_idx + 1, rec_list->count, x_vals->length);
        
        int counting = 0;
        int tau_current = 0;
        double d0_current = 0;
        int tau_count_this_pair = 0;
        
        for (int i = 0; i < x_vals->length; i++) {
            double dist = norm3d(x_vals->points[i], y_vals->points[i]);
            
            if (!counting) {
                if (dist >= soglia_min && dist <= soglia_max) {
                    counting = 1;
                    tau_current = 0;
                    d0_current = dist;
                }
            } else {
                tau_current++;
                if (dist >= delta_max) {
                    if (total_tau_count < 10000) {
                        all_tau_list[total_tau_count] = tau_current;
                        all_d0_list[total_tau_count] = d0_current;
                        total_tau_count++;
                        tau_count_this_pair++;
                    }
                    counting = 0;
                    tau_current = 0;
                }
            }
        }
        
        printf("      Trovati %d intervalli tau per questa ricorrenza\n", tau_count_this_pair);
    }
    
    // Calcola il Lyapunov finale usando tutti i tau raccolti
    double lyap_val;
    if (total_tau_count > 0) {
        double mean_tau = 0, mean_d0 = 0;
        for (int i = 0; i < total_tau_count; i++) {
            mean_tau += all_tau_list[i];
            mean_d0 += all_d0_list[i];
        }
        mean_tau /= total_tau_count;
        mean_d0 /= total_tau_count;
        
        lyap_val = (1.0 / (mean_tau * DT)) * log(delta_max / mean_d0);
        printf("  Totale intervalli tau trovati: %d\n", total_tau_count);
        printf("  Tau medio finale: %.2f, d0 medio finale: %.4f\n", mean_tau, mean_d0);
    } else {
        lyap_val = NAN;
        printf("  Nessun intervallo tau trovato in tutte le ricorrenze!\n");
    }
    
    free(all_tau_list);
    free(all_d0_list);
    return lyap_val;
}

// ===============================
// Funzione principale
// ===============================
int main() {
    // Inizializza il generatore di numeri casuali
    srand(time(NULL));
    
    printf("=== PARAMETRI SIMULAZIONE ===\n");
    printf("Numero di passi: %d\n", NUM_STEPS);
    printf("Passi minimi tra ricorrenze: %d\n", T_STEPS);
    printf("DELTA (soglia massima): %.1f\n", DELTA_MAX);
    printf("delta0 (soglia ricorrenze): %.1f\n", DELTA0);
    printf("Tolleranza: %.1f\n", TOL);
    printf("dt: %.2f\n", DT);
    printf("Metodo: %s\n", USE_EMBEDDING ? "Embedding 3D" : "Output diretto del modello");
    printf("Rumore: %s\n", ADD_NOISE ? "Attivo" : "Disattivo");
    if (ADD_NOISE) {
        printf("Intensità rumore: %.3f\n", NOISE_INTENSITY);
    }
    printf("\n");
    
    double *lyap_trad_list = malloc(NUM_RHO * sizeof(double));
    double *lyap_tau_list = malloc(NUM_RHO * sizeof(double));
    double *rho_values = malloc(NUM_RHO * sizeof(double));
    
    // Parametri embedding
    int embedding_tau = 10;
    int embedding_dim = 3;
    int component_to_embed = 0;  // 0=X, 1=Y, 2=Z
    
    if (USE_EMBEDDING) {
        printf("Parametri embedding: tau=%d, dim=%d, componente=%d\n", 
               embedding_tau, embedding_dim, component_to_embed);
    } else {
        printf("Usando output diretto del sistema di Lorenz (X, Y, Z)\n");
    }
    printf("\n");
    
    // Genera valori di rho (da 0 a 40)
    for (int i = 0; i < NUM_RHO; i++) {
        rho_values[i] = (double)i * 40.0 / (NUM_RHO - 1);
    }
    
    printf("=== INIZIO SIMULAZIONI ===\n");
    
    for (int i = 0; i < NUM_RHO; i++) {
        double rho = rho_values[i];
        printf("\n--- Simulazione %d/%d: rho = %.2f ---\n", i+1, NUM_RHO, rho);
        
        // Fissa il seed per riproducibilità
        srand(42 + i);
        
        // Genera traiettoria completa del sistema di Lorenz con rumore
        printf("Generando traiettoria %s del sistema di Lorenz...\n", 
               ADD_NOISE ? "con rumore" : "deterministica");
        Trajectory lorenz_traj;
        lorenz_traj.length = NUM_STEPS + 1;
        lorenz_traj.points = malloc(lorenz_traj.length * sizeof(Point3D));
        
        Point3D current = {1.0, 1.0, 1.0};
        lorenz_traj.points[0] = current;
        
        for (int step = 0; step < NUM_STEPS; step++) {
            if (step % 2000 == 0) {
                printf("  Passo %d/%d\n", step, NUM_STEPS);
            }
            current = lorenz_system_sde(current, rho, DT, NOISE_INTENSITY);
            lorenz_traj.points[step + 1] = current;
        }
        printf("Traiettoria completa generata: %d punti\n", lorenz_traj.length);
        
        // Scegli il metodo in base al flag
        Trajectory x_series;
        if (USE_EMBEDDING) {
            printf("METODO EMBEDDING:\n");
            // Estrai la componente desiderata
            double *component_series = malloc(lorenz_traj.length * sizeof(double));
            for (int j = 0; j < lorenz_traj.length; j++) {
                switch(component_to_embed) {
                    case 0: component_series[j] = lorenz_traj.points[j].x; break;
                    case 1: component_series[j] = lorenz_traj.points[j].y; break;
                    case 2: component_series[j] = lorenz_traj.points[j].z; break;
                }
            }
            printf("Componente %s estratta: %d punti\n", 
                   (component_to_embed == 0) ? "X" : (component_to_embed == 1) ? "Y" : "Z", 
                   lorenz_traj.length);
            
            // Crea embedding 3D
            x_series = create_3d_embedding(component_series, lorenz_traj.length, 
                                         embedding_tau, embedding_dim);
            free(component_series);
            
            if (x_series.length <= 0) {
                lyap_trad_list[i] = NAN;
                lyap_tau_list[i] = NAN;
                free_trajectory(&lorenz_traj);
                continue;
            }
        } else {
            printf("METODO DIRETTO:\n");
            x_series = lorenz_traj;  // Usa direttamente la traiettoria 3D
            printf("Usando traiettoria diretta 3D: %d punti\n", x_series.length);
        }
        
        // Lyapunov tradizionale
        printf("Calcolando Lyapunov tradizionale (sistema deterministico equivalente)...\n");
        Point3D v = {1.0, 0.0, 0.0};
        double lyap_sum = 0;
        
        for (int j = 0; j < lorenz_traj.length; j++) {
            Matrix *J = jacobian(lorenz_traj.points[j], rho);
            Point3D Jv = matrix_vector_mult(J, v);
            v = add_points(v, scale_point(Jv, DT));
            double norm_v = vector_norm(v);
            v = normalize_vector(v);
            lyap_sum += log(norm_v);
            free_matrix(J);
        }
        
        double lyap_trad = lyap_sum / (NUM_STEPS * DT);
        lyap_trad_list[i] = lyap_trad;
        printf("Lyapunov tradizionale: %.4f\n", lyap_trad);
        
        // Lyap_tau
        printf("Calcolando Lyapunov tau su metodo %s %s...\n", 
               USE_EMBEDDING ? "embedding" : "diretto",
               ADD_NOISE ? "con rumore" : "senza rumore");
        RecurrenceList recurrence_pairs = select_y_by_recurrence(&x_series, DELTA0, TOL, T_STEPS, MAX_RECURRENCES);
        
        if (recurrence_pairs.count > 0) {
            double lyap_tau_val = lyap_tau_continuous_multi(&recurrence_pairs, DELTA_MAX, DELTA0, TOL);
            lyap_tau_list[i] = lyap_tau_val;
        } else {
            lyap_tau_list[i] = NAN;
            printf("  Lyapunov tau: Non calcolabile (nessuna ricorrenza)\n");
        }
        
        printf("RISULTATO: rho=%.2f | Lyap trad=%.3f | Lyap tau=%.3f\n", 
               rho, lyap_trad, lyap_tau_list[i]);
        
        // Cleanup per questa iterazione
        if (USE_EMBEDDING) {
            free_trajectory(&x_series);
        }
        free_trajectory(&lorenz_traj);
        free_recurrence_list(&recurrence_pairs);
    }
    
    // Salva risultati in file
    printf("\n=== SALVATAGGIO RISULTATI ===\n");
    const char *method_str = USE_EMBEDDING ? "embedding" : "direct";
    const char *noise_str = ADD_NOISE ? "noise" : "no_noise";
    char filename[256];
    snprintf(filename, sizeof(filename), "lyapunov_vs_rho_%s_%s.dat", method_str, noise_str);
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "# rho lyap_trad lyap_tau\n");
    for (int i = 0; i < NUM_RHO; i++) {
        fprintf(fp, "%.6f %.6f %.6f\n", rho_values[i], lyap_trad_list[i], lyap_tau_list[i]);
    }
    fclose(fp);
    printf("Risultati salvati in '%s'\n", filename);
    
    // Cleanup finale
    free(lyap_trad_list);
    free(lyap_tau_list);
    free(rho_values);
    
    printf("=== SIMULAZIONE COMPLETATA ===\n");
    return 0;
}
