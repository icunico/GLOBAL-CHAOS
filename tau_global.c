#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>

// ===== NetCDF error handling macro =====
#define NC_CHECK(stat) { if((stat) != NC_NOERR) { fprintf(stderr, "NetCDF error: %s at line %d\n", nc_strerror(stat), __LINE__); exit(1); } }

// ===== Simulation parameters =====
#define T_STEPS 90          // Tempo tra ricorrenze
#define DELTAMAX 20.0       // Soglia massima
#define DT 1.0              // Time step (daily data)
#define TOL 0.4              // Tolleranza
#define DELTA0 0.05          // Delta0 = 0.05

// ===== Flags =====
#define ADD_NOISE 0
#define ZERO_DETERMINISTIC 0
#define USE_EMBEDDING 1
#define EMB_DIM 3
#define TAU_EMB 10

// ===== Structs =====
typedef struct { 
    double x[EMB_DIM]; 
} VecEmb;

typedef struct {
    VecEmb *x_vals;
    VecEmb *y_vals;
    int length;
} RecurrencePair;

// ===== Funzioni utili =====
double randn() {
    double u1 = ((double)rand()) / RAND_MAX;
    double u2 = ((double)rand()) / RAND_MAX;
    return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

double norm(VecEmb a, VecEmb b) {
    double s = 0;
    for(int i=0; i<EMB_DIM; i++){
        double d = a.x[i] - b.x[i];
        s += d*d;
    }
    return sqrt(s);
}

// ===== Create embedding vector =====
VecEmb create_embedding(double *series, int start_idx, int len) {
    VecEmb v;
    for(int k = 0; k < EMB_DIM; k++) {
        int idx = start_idx + k * TAU_EMB;
        if(idx < len) {
            v.x[k] = series[idx];
        } else {
            v.x[k] = series[len - 1];
        }
    }
    return v;
}

// ===== Ricorrenze con embedding =====
RecurrencePair* select_y_by_recurrence(double *series, int N, int *n_pairs, 
                                       double delta0) {
    double soglia_min = delta0 * (1 - TOL);
    double soglia_max = delta0 * (1 + TOL);
    int max_rec = 50;
    RecurrencePair *pairs = malloc(max_rec * sizeof(RecurrencePair));
    int rec_found = 0;
    int last_t1_used = -T_STEPS;

    int max_index = N - (USE_EMBEDDING ? (EMB_DIM - 1) * TAU_EMB : 0);
    if(max_index <= 0) {
        *n_pairs = 0;
        return pairs;
    }

    for(int t1 = 0; t1 < max_index - T_STEPS; t1++) {
        if(t1 - last_t1_used < T_STEPS) continue;
        
        for(int t2 = t1 + T_STEPS; t2 < max_index; t2++) {
            VecEmb v1 = create_embedding(series, t1, N);
            VecEmb v2 = create_embedding(series, t2, N);

            double d = norm(v1, v2);
            if(d >= soglia_min && d <= soglia_max) {
                int max_len = N - (t1 > t2 ? t1 : t2) - (EMB_DIM - 1) * TAU_EMB;
                if(max_len <= 0) continue;
                
                pairs[rec_found].x_vals = malloc(max_len * sizeof(VecEmb));
                pairs[rec_found].y_vals = malloc(max_len * sizeof(VecEmb));
                
                for(int k = 0; k < max_len; k++) {
                    pairs[rec_found].x_vals[k] = create_embedding(series, t1 + k, N);
                    pairs[rec_found].y_vals[k] = create_embedding(series, t2 + k, N);
                }
                
                pairs[rec_found].length = max_len;
                
                rec_found++;
                last_t1_used = t1;
                if(rec_found >= max_rec) break;
                break;
            }
        }
        if(rec_found >= max_rec) break;
    }
    
    *n_pairs = rec_found;
    return pairs;
}

// ===== Calcolo Lyapunov tau (Palladin method) =====
double lyap_palladin(double *series, int N, double delta0) {
    int n_pairs;
    RecurrencePair *pairs = select_y_by_recurrence(series, N, &n_pairs, delta0);
    
    if(n_pairs == 0) {
        for(int p = 0; p < n_pairs; p++) {
            free(pairs[p].x_vals);
            free(pairs[p].y_vals);
        }
        free(pairs);
        return NAN;
    }
    
    int tau_count = 0;
    double tau_sum = 0, d0_sum = 0;
    double soglia_min = delta0 * (1 - TOL);
    double soglia_max = delta0 * (1 + TOL);

    for(int p = 0; p < n_pairs; p++) {
        int counting = 0, tau_current = 0;
        double d0_current = 0;
        
        for(int i = 0; i < pairs[p].length; i++) {
            double d = norm(pairs[p].x_vals[i], pairs[p].y_vals[i]);
            
            if(!counting) {
                if(d >= soglia_min && d <= soglia_max) {
                    counting = 1; 
                    tau_current = 0; 
                    d0_current = d;
                }
            } else {
                tau_current++;
                if(d >= DELTAMAX) {
                    tau_sum += tau_current;
                    d0_sum += d0_current;
                    tau_count++;
                    counting = 0;
                }
            }
        }
    }

    for(int p = 0; p < n_pairs; p++) {
        free(pairs[p].x_vals);
        free(pairs[p].y_vals);
    }
    free(pairs);

    if(tau_count == 0) return NAN;

    double mean_tau = tau_sum / tau_count;
    double mean_d0 = d0_sum / tau_count;
    double lyap = log(DELTAMAX / mean_d0) / (mean_tau * DT);
    
    return lyap;
}

// ===== Wolf Lyapunov =====
double wolf_lyapunov(double *series, int N) {
    if(N < (EMB_DIM - 1) * TAU_EMB + 10) return NAN;
    
    double evolution_time = 20;
    double sum_lyap = 0;
    int count = 0;
    
    for(int start = 0; start < N - evolution_time - (EMB_DIM-1)*TAU_EMB; 
        start += evolution_time/2) {
        
        VecEmb ref = create_embedding(series, start, N);
        
        double min_dist = INFINITY;
        int nearest_idx = -1;
        
        for(int i = 0; i < N - (EMB_DIM-1)*TAU_EMB; i++) {
            if(abs(i - start) < T_STEPS) continue;
            
            VecEmb test = create_embedding(series, i, N);
            double dist = norm(ref, test);
            
            if(dist < min_dist && dist > 0) {
                min_dist = dist;
                nearest_idx = i;
            }
        }
        
        if(nearest_idx >= 0 && min_dist > 0) {
            VecEmb ref_evolved = create_embedding(series, start + evolution_time, N);
            VecEmb neigh_evolved = create_embedding(series, nearest_idx + evolution_time, N);
            
            double final_dist = norm(ref_evolved, neigh_evolved);
            
            if(final_dist > 0 && min_dist > 0) {
                sum_lyap += log(final_dist / min_dist) / evolution_time;
                count++;
            }
        }
    }
    
    return (count > 0) ? sum_lyap / count : NAN;
}

// ===== MAIN =====
int main(int argc, char **argv) {
    if(argc < 2) {
        printf("Usage: %s <chlorophyll_netcdf_file>\n", argv[0]);
        return 1;
    }
    
    // Apri file NetCDF
    int ncid;
    NC_CHECK(nc_open(argv[1], NC_NOWRITE, &ncid));
    
    // Ottieni dimensioni
    int dimid_time, dimid_lat, dimid_lon;
    size_t n_time, n_lat, n_lon;
    
    NC_CHECK(nc_inq_dimid(ncid, "time", &dimid_time));
    NC_CHECK(nc_inq_dimid(ncid, "latitude", &dimid_lat));
    NC_CHECK(nc_inq_dimid(ncid, "longitude", &dimid_lon));
    
    NC_CHECK(nc_inq_dimlen(ncid, dimid_time, &n_time));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_lat, &n_lat));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_lon, &n_lon));
    
    // Ottieni variabile CHL
    int varid_chl;
    NC_CHECK(nc_inq_varid(ncid, "CHL", &varid_chl));
    
    // Leggi coordinate
    float *latitudes = malloc(n_lat * sizeof(float));
    float *longitudes = malloc(n_lon * sizeof(float));
    
    int varid_lat, varid_lon;
    NC_CHECK(nc_inq_varid(ncid, "latitude", &varid_lat));
    NC_CHECK(nc_inq_varid(ncid, "longitude", &varid_lon));
    
    NC_CHECK(nc_get_var_float(ncid, varid_lat, latitudes));
    NC_CHECK(nc_get_var_float(ncid, varid_lon, longitudes));
    
    printf("=== INIZIO SIMULAZIONE ===\n");
    printf("File: %s\n", argv[1]);
    printf("Dimensioni GLOBALI: time=%zu, lat=%zu, lon=%zu\n", n_time, n_lat, n_lon);
    printf("Parametri: delta0=%.3f, T_STEPS=%d, EMB_DIM=%d, TAU_EMB=%d\n", 
           DELTA0, T_STEPS, EMB_DIM, TAU_EMB);
    printf("Processo TUTTA la griglia globale (%zu x %zu = %zu punti)\n", 
           n_lat, n_lon, n_lat * n_lon);
    
    // USA TUTTA LA GRIGLIA (nessuna selezione)
    int nlat_sel = n_lat;
    int nlon_sel = n_lon;
    int npoints = nlat_sel * nlon_sel;
    
    printf("Totale punti da processare: %d\n", npoints);
    printf("ATTENZIONE: Potrebbe richiedere molto tempo e molta memoria!\n");
    
    // Alloca array per risultati
    double **lyap_wolf_all = malloc(nlat_sel * sizeof(double*));
    double **lyap_palladin_all = malloc(nlat_sel * sizeof(double*));
    
    for(int i = 0; i < nlat_sel; i++) {
        lyap_wolf_all[i] = malloc(nlon_sel * sizeof(double));
        lyap_palladin_all[i] = malloc(nlon_sel * sizeof(double));
    }
    
    // Buffer per serie temporale
    double *series_buffer = malloc(n_time * sizeof(double));
    
    // Loop su TUTTI i punti della griglia
    int point_count = 0;
    time_t start_time = time(NULL);
    
    for(int ilat = 0; ilat < nlat_sel; ilat++) {
        for(int ilon = 0; ilon < nlon_sel; ilon++) {
            
            point_count++;
            
            // Stima tempo rimanente ogni 1000 punti
            if(point_count % 1000 == 0 || point_count == 1) {
                time_t current_time = time(NULL);
                double elapsed = difftime(current_time, start_time);
                double points_per_sec = point_count / elapsed;
                int remaining_points = npoints - point_count;
                double remaining_hours = (remaining_points / points_per_sec) / 3600;
                
                printf("\rProcessando punto %d/%d (lat=%d, lon=%d) - %.1f punti/sec - %.1f ore rimaste", 
                       point_count, npoints, ilat, ilon, points_per_sec, remaining_hours);
                fflush(stdout);
            }
            
            // Leggi serie temporale per questo punto
            size_t start[3] = {0, (size_t)ilat, (size_t)ilon};
            size_t count[3] = {n_time, 1, 1};
            
            NC_CHECK(nc_get_vara_double(ncid, varid_chl, start, count, series_buffer));
            
            // Filtra dati validi
            int valid_count = 0;
            for(size_t t = 0; t < n_time; t++) {
                if(!isnan(series_buffer[t]) && series_buffer[t] > 0 && series_buffer[t] < 100) {
                    series_buffer[valid_count++] = series_buffer[t];
                }
            }
            
            if(valid_count < 100) {
                lyap_wolf_all[ilat][ilon] = NAN;
                lyap_palladin_all[ilat][ilon] = NAN;
                continue;
            }
            
            // Calcola Wolf Lyapunov
            lyap_wolf_all[ilat][ilon] = wolf_lyapunov(series_buffer, valid_count);
            
            // Calcola Palladin con delta0=0.05
            lyap_palladin_all[ilat][ilon] = lyap_palladin(series_buffer, valid_count, DELTA0);
        }
    }
    printf("\n\n");
    
    // Chiudi file NetCDF
    NC_CHECK(nc_close(ncid));
    
    // Salva risultati
    
    // 1. Coordinate
    FILE *f_coords = fopen("coordinate_globali.dat", "w");
    fprintf(f_coords, "# lat lon\n");
    for(int ilat = 0; ilat < nlat_sel; ilat++) {
        for(int ilon = 0; ilon < nlon_sel; ilon++) {
            fprintf(f_coords, "%.4f %.4f\n", latitudes[ilat], longitudes[ilon]);
        }
    }
    fclose(f_coords);
    
    // 2. Wolf Lyapunov
    FILE *f_wolf = fopen("lyapunov_wolf_globale.dat", "w");
    fprintf(f_wolf, "# Wolf Lyapunov - griglia globale\n");
    for(int ilat = 0; ilat < nlat_sel; ilat++) {
        for(int ilon = 0; ilon < nlon_sel; ilon++) {
            fprintf(f_wolf, "%.6f\n", lyap_wolf_all[ilat][ilon]);
        }
    }
    fclose(f_wolf);
    
    // 3. Palladin Lyapunov (delta0=0.05)
    FILE *f_pall = fopen("lyapunov_palladin_globale.dat", "w");
    fprintf(f_pall, "# Palladin Lyapunov with delta0=%.3f - griglia globale\n", DELTA0);
    for(int ilat = 0; ilat < nlat_sel; ilat++) {
        for(int ilon = 0; ilon < nlon_sel; ilon++) {
            fprintf(f_pall, "%.6f\n", lyap_palladin_all[ilat][ilon]);
        }
    }
    fclose(f_pall);
    
    // 4. Metadata
    FILE *f_meta = fopen("metadata_globali.txt", "w");
    fprintf(f_meta, "DELTA0=%.3f\n", DELTA0);
    fprintf(f_meta, "T_STEPS=%d\n", T_STEPS);
    fprintf(f_meta, "EMB_DIM=%d\n", EMB_DIM);
    fprintf(f_meta, "TAU_EMB=%d\n", TAU_EMB);
    fprintf(f_meta, "NLAT=%d\n", nlat_sel);
    fprintf(f_meta, "NLON=%d\n", nlon_sel);
    fprintf(f_meta, "N_TIME=%zu\n", n_time);
    fprintf(f_meta, "TOT_POINTS=%d\n", npoints);
    fclose(f_meta);
    
    // Libera memoria
    free(series_buffer);
    free(latitudes);
    free(longitudes);
    
    for(int i = 0; i < nlat_sel; i++) {
        free(lyap_wolf_all[i]);
        free(lyap_palladin_all[i]);
    }
    free(lyap_wolf_all);
    free(lyap_palladin_all);
    
    printf("=== SIMULAZIONE GLOBALE COMPLETATA ===\n");
    printf("File generati:\n");
    printf("  - coordinate_globali.dat\n");
    printf("  - lyapunov_wolf_globale.dat\n");
    printf("  - lyapunov_palladin_globale.dat\n");
    printf("  - metadata_globali.txt\n");
    printf("\nDelta0 utilizzato: %.3f\n", DELTA0);
    printf("Tempo totale: %.1f secondi\n", difftime(time(NULL), start_time));
    
    return 0;
}
