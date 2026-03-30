#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <sys/time.h>

// ===== Parametri =====
#define DELTAMAX 0.2
#define DELTA0 0.02
#define DT 0.01
#define TOL 0.4
#define T_STEPS 90
#define EMB_DIM 3
#define TAU_EMB 10

// ===== Struct per vettori embeddati =====
typedef struct {
    double x[EMB_DIM];
} VecEmb;

// ===== Funzione norma =====
double norm_vec(const VecEmb* a, const VecEmb* b) {
    double s = 0;
    for(int i = 0; i < EMB_DIM; i++) {
        double d = a->x[i] - b->x[i];
        s += d * d;
    }
    return sqrt(s);
}

// ===== Creazione embedding =====
VecEmb* create_embedding(const double* series, int N, int* emb_length) {
    int max_index = N - (EMB_DIM - 1) * TAU_EMB;
    if(max_index <= 0) {
        *emb_length = 0;
        return NULL;
    }
    
    VecEmb* embedded = (VecEmb*)malloc(max_index * sizeof(VecEmb));
    for(int i = 0; i < max_index; i++) {
        for(int dim = 0; dim < EMB_DIM; dim++) {
            int idx = i + dim * TAU_EMB;
            embedded[i].x[dim] = series[idx];
        }
    }
    
    *emb_length = max_index;
    return embedded;
}

// ===== Lyapunov di Paladin (SOLO QUESTO) =====
double lyapunov_paladin(const double* series, int N) {
    if(N < 100) return NAN;
    
    int N_emb;
    VecEmb* series_emb = create_embedding(series, N, &N_emb);
    if(!series_emb || N_emb < T_STEPS + 10) {
        free(series_emb);
        return NAN;
    }
    
    double soglia_min = DELTA0 * (1 - TOL);
    double soglia_max = DELTA0 * (1 + TOL);
    
    // Trova ricorrenze
    int max_pairs = 50;
    int n_pairs = 0;
    int* t1_list = (int*)malloc(max_pairs * sizeof(int));
    int* t2_list = (int*)malloc(max_pairs * sizeof(int));
    int last_t1_used = -T_STEPS;
    
    for(int t1 = 0; t1 < N_emb - T_STEPS; t1++) {
        if(t1 - last_t1_used < T_STEPS) continue;
        
        for(int t2 = t1 + T_STEPS; t2 < N_emb; t2++) {
            double d = norm_vec(&series_emb[t1], &series_emb[t2]);
            
            if(d >= soglia_min && d <= soglia_max) {
                t1_list[n_pairs] = t1;
                t2_list[n_pairs] = t2;
                n_pairs++;
                last_t1_used = t1;
                break;
            }
        }
        if(n_pairs >= max_pairs) break;
    }
    
    if(n_pairs == 0) {
        free(series_emb);
        free(t1_list);
        free(t2_list);
        return NAN;
    }
    
    // Calcola tempo di predicibilità
    int tau_count = 0;
    double tau_sum = 0, d0_sum = 0;
    
    for(int p = 0; p < n_pairs; p++) {
        int t1 = t1_list[p];
        int t2 = t2_list[p];
        int max_len = N_emb - (t1 > t2 ? t1 : t2);
        
        int counting = 0, tau_current = 0;
        double d0_current = 0;
        
        for(int i = 0; i < max_len; i++) {
            double d = norm_vec(&series_emb[t1 + i], &series_emb[t2 + i]);
            
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
    
    free(series_emb);
    free(t1_list);
    free(t2_list);
    
    if(tau_count == 0) return NAN;
    
    double mean_tau = tau_sum / tau_count;
    double mean_d0 = d0_sum / tau_count;
    double Tcrit = mean_tau * DT;
    double Tpred = Tcrit * log(DELTAMAX / DELTA0) / log(DELTAMAX / mean_d0);
    
    return Tpred;
}

// ===== Pulizia serie da NaN =====
int clean_series(const float* input, int N, double** output) {
    int count = 0;
    for(int i = 0; i < N; i++) {
        if(!isnan(input[i]) && input[i] > 0 && input[i] < 100) {
            count++;
        }
    }
    
    if(count < 100) return 0;
    
    *output = (double*)malloc(count * sizeof(double));
    int idx = 0;
    for(int i = 0; i < N; i++) {
        if(!isnan(input[i]) && input[i] > 0 && input[i] < 100) {
            (*output)[idx++] = (double)input[i];
        }
    }
    
    return count;
}

// ===== Leggi informazioni dal netCDF =====
int get_netcdf_info(const char* filename, const char* varname,
                    int* nlat, int* nlon, int* ntime,
                    double** lats, double** lons) {
    
    int ncid;
    int retval;
    
    if ((retval = nc_open(filename, NC_NOWRITE, &ncid))) {
        printf("ERRORE: Impossibile aprire %s: %s\n", filename, nc_strerror(retval));
        return 0;
    }
    
    int lat_dimid, lon_dimid, time_dimid;
    if (nc_inq_dimid(ncid, "latitude", &lat_dimid) != NC_NOERR ||
        nc_inq_dimid(ncid, "longitude", &lon_dimid) != NC_NOERR ||
        nc_inq_dimid(ncid, "time", &time_dimid) != NC_NOERR) {
        printf("ERRORE: Dimensioni non trovate\n");
        nc_close(ncid);
        return 0;
    }
    
    nc_inq_dimlen(ncid, lat_dimid, (size_t*)nlat);
    nc_inq_dimlen(ncid, lon_dimid, (size_t*)nlon);
    nc_inq_dimlen(ncid, time_dimid, (size_t*)ntime);
    
    int lat_varid, lon_varid;
    nc_inq_varid(ncid, "latitude", &lat_varid);
    nc_inq_varid(ncid, "longitude", &lon_varid);
    
    *lats = (double*)malloc(*nlat * sizeof(double));
    *lons = (double*)malloc(*nlon * sizeof(double));
    
    nc_get_var_double(ncid, lat_varid, *lats);
    nc_get_var_double(ncid, lon_varid, *lons);
    
    nc_close(ncid);
    return 1;
}

// ===== Leggi un punto specifico =====
float* read_netcdf_point(const char* filename, const char* varname,
                         int lat_idx, int lon_idx, int ntime) {
    
    int ncid;
    nc_open(filename, NC_NOWRITE, &ncid);
    
    int varid;
    nc_inq_varid(ncid, varname, &varid);
    
    float* data = (float*)malloc(ntime * sizeof(float));
    size_t start[3] = {0, (size_t)lat_idx, (size_t)lon_idx};
    size_t count[3] = {(size_t)ntime, 1, 1};
    
    nc_get_vara_float(ncid, varid, start, count, data);
    nc_close(ncid);
    
    return data;
}

// ===== Funzione per il tempo =====
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// ===== MAIN =====
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 3) {
        if (rank == 0) {
            printf("Uso: mpirun -np N %s <file.nc> <varname>\n", argv[0]);
            printf("Calcola SOLO Paladin (tempo di predicibilità) su tutta l'area\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    const char* filename = argv[1];
    const char* varname = argv[2];
    
    int nlat, nlon, ntime;
    double *lats = NULL, *lons = NULL;
    
    double t_start_total = get_time();
    
    // Solo rank 0 legge info
    if (rank == 0) {
        printf("\n=== Lyapunov Paladin - AREA COMPLETA ===\n");
        printf("File: %s\n", filename);
        printf("Variabile: %s\n", varname);
        printf("Processi MPI: %d\n", size);
        printf("Calcolo: SOLO Paladin (NO Wolf)\n");
        
        if (!get_netcdf_info(filename, varname, &nlat, &nlon, &ntime, &lats, &lons)) {
            printf("ERRORE: Impossibile leggere il file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("Dimensioni: lat=%d, lon=%d, time=%d\n", nlat, nlon, ntime);
        printf("Totale punti: %d\n", nlat * nlon);
        printf("================================\n\n");
    }
    
    // Broadcast dimensioni
    MPI_Bcast(&nlat, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nlon, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        lats = (double*)malloc(nlat * sizeof(double));
        lons = (double*)malloc(nlon * sizeof(double));
    }
    
    MPI_Bcast(lats, nlat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lons, nlon, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Distribuisci TUTTI i punti dell'area
    int total_points = nlat * nlon;
    int points_per_rank = total_points / size;
    int remainder = total_points % size;
    
    int start_idx, local_count;
    if (rank < remainder) {
        start_idx = rank * (points_per_rank + 1);
        local_count = points_per_rank + 1;
    } else {
        start_idx = rank * points_per_rank + remainder;
        local_count = points_per_rank;
    }
    
    printf("[Rank %d] Processando %d punti su %d totali (dal %d al %d)\n", 
           rank, local_count, total_points, start_idx, start_idx + local_count - 1);
    
    // Alloca risultati locali (SOLO Paladin)
    double* local_paladin = (double*)malloc(local_count * sizeof(double));
    double* local_lats = (double*)malloc(local_count * sizeof(double));
    double* local_lons = (double*)malloc(local_count * sizeof(double));
    
    double t_start_calc = get_time();
    int local_valid = 0;
    
    // Processa ogni punto
    for (int i = 0; i < local_count; i++) {
        int global_idx = start_idx + i;
        int lat_idx = global_idx / nlon;
        int lon_idx = global_idx % nlon;
        
        local_lats[i] = lats[lat_idx];
        local_lons[i] = lons[lon_idx];
        
        // Leggi serie
        float* raw_data = read_netcdf_point(filename, varname, lat_idx, lon_idx, ntime);
        
        // Pulisci
        double* clean_data;
        int clean_len = clean_series(raw_data, ntime, &clean_data);
        free(raw_data);
        
        if (clean_len >= 100) {
            local_paladin[i] = lyapunov_paladin(clean_data, clean_len);
            free(clean_data);
            if (!isnan(local_paladin[i])) local_valid++;
        } else {
            local_paladin[i] = NAN;
        }
        
        // Progress report ogni 100 punti
        if (rank == 0 && (i+1) % 100 == 0) {
            double t_now = get_time();
            double elapsed = t_now - t_start_calc;
            double est_totale = (elapsed / (i+1)) * local_count;
            printf("[Rank 0] Progresso: %d/%d punti (%.1f%%) - tempo stimato: %.1f min\n", 
                   i+1, local_count, 100.0*(i+1)/local_count, est_totale/60.0);
        }
    }
    
    double t_end_calc = get_time();
    
    printf("[Rank %d] Completato: %d/%d punti validi\n", 
           rank, local_valid, local_count);
    
    // Prepara buffer per gather
    int* recvcounts = NULL;
    int* displs = NULL;
    double* global_paladin = NULL;
    double* global_lats_all = NULL;
    double* global_lons_all = NULL;
    
    if (rank == 0) {
        recvcounts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));
        
        int offset = 0;
        for (int r = 0; r < size; r++) {
            int r_count = (r < remainder) ? points_per_rank + 1 : points_per_rank;
            recvcounts[r] = r_count;
            displs[r] = offset;
            offset += r_count;
        }
        
        global_paladin = (double*)malloc(total_points * sizeof(double));
        global_lats_all = (double*)malloc(total_points * sizeof(double));
        global_lons_all = (double*)malloc(total_points * sizeof(double));
    }
    
    // Raccogli risultati
    MPI_Gatherv(local_paladin, local_count, MPI_DOUBLE,
                global_paladin, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(local_lats, local_count, MPI_DOUBLE,
                global_lats_all, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(local_lons, local_count, MPI_DOUBLE,
                global_lons_all, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    double t_end_total = get_time();
    
    // Rank 0 salva risultati
    if (rank == 0) {
        printf("\n=== RISULTATI AREA COMPLETA ===\n");
        printf("Tempo totale: %.2f secondi (%.1f minuti)\n", 
               t_end_total - t_start_total, (t_end_total - t_start_total)/60.0);
        printf("Tempo calcolo: %.2f secondi\n", t_end_calc - t_start_calc);
        
        int valid_count = 0;
        for (int i = 0; i < total_points; i++) {
            if (!isnan(global_paladin[i])) valid_count++;
        }
        
        printf("Punti validi: %d/%d (%.1f%%)\n", 
               valid_count, total_points, 100.0 * valid_count / total_points);
        
        // Salva in formato DAT
        FILE* f = fopen("lyapunov_paladin_area.dat", "w");
        fprintf(f, "# latitude\tlongitude\tpaladin_tau\n");
        for (int i = 0; i < total_points; i++) {
            fprintf(f, "%f\t%f\t%f\n", 
                    global_lats_all[i], global_lons_all[i], global_paladin[i]);
        }
        fclose(f);
        
        printf("\nRisultati salvati in: lyapunov_paladin_area.dat\n");
        
        // Statistiche
        double sum = 0;
        int n = 0;
        double min_val = 1e30, max_val = -1e30;
        for (int i = 0; i < total_points; i++) {
            if (!isnan(global_paladin[i])) {
                sum += global_paladin[i];
                n++;
                if (global_paladin[i] < min_val) min_val = global_paladin[i];
                if (global_paladin[i] > max_val) max_val = global_paladin[i];
            }
        }
        if (n > 0) {
            printf("\nStatistiche τ (tempo di predicibilità):\n");
            printf("  Media: %.6f\n", sum/n);
            printf("  Min: %.6f\n", min_val);
            printf("  Max: %.6f\n", max_val);
        }
        
        free(recvcounts);
        free(displs);
        free(global_paladin);
        free(global_lats_all);
        free(global_lons_all);
    }
    
    // Libera memoria
    free(local_paladin);
    free(local_lats);
    free(local_lons);
    free(lats);
    free(lons);
    
    MPI_Finalize();
    return 0;
}
