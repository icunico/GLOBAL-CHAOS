#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>

// ===== Parametri di default =====
#define DELTAMAX_DEFAULT 0.5  // Valore di fallback se CHL_lognormal non disponibile
#define DELTA0 0.02
#define DT 1.0
#define TOL 0.4
#define T_STEPS 90
#define EMB_DIM 3
#define TAU_EMB 10
#define EVOLVE_STEPS 20000

// ===== DECIMAZIONE SPAZIALE =====
#define SPATIAL_DECIMATION 2  // Prendi 1 punto ogni 5 (sia lat che lon)

// ===== Struct per vettori embeddati =====
typedef struct {
    double x[EMB_DIM];
    int valid;  // Flag per indicare se il vettore è valido (senza NaN)
} VecEmb;

// ===== Struct per griglia spaziale dei parametri =====
typedef struct {
    int nlat, nlon, ntime;
    double *lats, *lons;
    float *deltamax_grid;  // Griglia spaziale di DELTAMAX
    int valid_deltamax;    // Flag se la griglia è stata caricata con successo
} SpatialParams;

// ===== Funzione norma (solo per vettori validi) =====
double norm_vec(const VecEmb* a, const VecEmb* b) {
    double s = 0;
    for(int i = 0; i < EMB_DIM; i++) {
        double d = a->x[i] - b->x[i];
        s += d * d;
    }
    return sqrt(s);
}

// ===== Verifica se un vettore contiene NaN =====
int has_nan(const VecEmb* v) {
    for(int i = 0; i < EMB_DIM; i++) {
        if(isnan(v->x[i])) return 1;
    }
    return 0;
}

// ===== Carica griglia DELTAMAX da NetCDF =====
int load_deltamax_grid(const char* deltamax_filename, SpatialParams* params) {
    int ncid, retval;
    
    printf("Tentativo di caricamento DELTAMAX da: %s\n", deltamax_filename);
    
    if ((retval = nc_open(deltamax_filename, NC_NOWRITE, &ncid))) {
        printf("⚠ ATTENZIONE: Impossibile aprire %s: %s\n", 
               deltamax_filename, nc_strerror(retval));
        printf("   Userò DELTAMAX costante = %f\n", DELTAMAX_DEFAULT);
        params->valid_deltamax = 0;
        return 0;
    }
    
    // Verifica se esiste la variabile CHL_lognormal
    int varid;
    if ((retval = nc_inq_varid(ncid, "CHL_lognormal", &varid))) {
        printf("⚠ ATTENZIONE: Variabile 'CHL_lognormal' non trovata in %s: %s\n", 
               deltamax_filename, nc_strerror(retval));
        printf("   Userò DELTAMAX costante = %f\n", DELTAMAX_DEFAULT);
        nc_close(ncid);
        params->valid_deltamax = 0;
        return 0;
    }
    
    // Ottieni dimensioni (dovrebbero essere compatibili con il file di input)
    int lat_dimid, lon_dimid;
    if (nc_inq_dimid(ncid, "latitude", &lat_dimid) != NC_NOERR ||
        nc_inq_dimid(ncid, "longitude", &lon_dimid) != NC_NOERR) {
        printf("⚠ ERRORE: Dimensioni latitude/longitude non trovate in %s\n", deltamax_filename);
        nc_close(ncid);
        params->valid_deltamax = 0;
        return 0;
    }
    
    size_t nlat_deltamax, nlon_deltamax;
    nc_inq_dimlen(ncid, lat_dimid, &nlat_deltamax);
    nc_inq_dimlen(ncid, lon_dimid, &nlon_deltamax);
    
    // Verifica compatibilità dimensioni
    if (nlat_deltamax != params->nlat || nlon_deltamax != params->nlon) {
        printf("⚠ ERRORE: Dimensioni incompatibili!\n");
        printf("   Input: lat=%d, lon=%d\n", params->nlat, params->nlon);
        printf("   DELTAMAX: lat=%zu, lon=%zu\n", nlat_deltamax, nlon_deltamax);
        nc_close(ncid);
        params->valid_deltamax = 0;
        return 0;
    }
    
    // Alloca e leggi i dati
    params->deltamax_grid = (float*)malloc(params->nlat * params->nlon * sizeof(float));
    
    if ((retval = nc_get_var_float(ncid, varid, params->deltamax_grid))) {
        printf("⚠ ERRORE: Impossibile leggere CHL_lognormal: %s\n", nc_strerror(retval));
        free(params->deltamax_grid);
        params->deltamax_grid = NULL;
        nc_close(ncid);
        params->valid_deltamax = 0;
        return 0;
    }
    
    nc_close(ncid);
    
    // Analizza i valori caricati
    int valid_count = 0, nan_count = 0, negative_count = 0;
    float min_val = 1e30, max_val = -1e30;
    double sum = 0;
    
    for (int i = 0; i < params->nlat * params->nlon; i++) {
        if (isnan(params->deltamax_grid[i])) {
            nan_count++;
        } else if (params->deltamax_grid[i] <= 0) {
            negative_count++;
        } else {
            valid_count++;
            sum += params->deltamax_grid[i];
            if (params->deltamax_grid[i] < min_val) min_val = params->deltamax_grid[i];
            if (params->deltamax_grid[i] > max_val) max_val = params->deltamax_grid[i];
        }
    }
    
    printf("✓ DELTAMAX caricato con successo da CHL_lognormal:\n");
    printf("   Punti validi: %d/%d (%.1f%%)\n", 
           valid_count, params->nlat * params->nlon, 
           100.0 * valid_count / (params->nlat * params->nlon));
    printf("   NaN: %d, Negativi/Zero: %d\n", nan_count, negative_count);
    
    if (valid_count > 0) {
        printf("   Range: [%.6f, %.6f]\n", min_val, max_val);
        printf("   Media: %.6f\n", sum / valid_count);
    }
    
    params->valid_deltamax = 1;
    return 1;
}

// ===== Ottieni DELTAMAX per un punto specifico =====
double get_deltamax_for_point(const SpatialParams* params, int lat_idx, int lon_idx) {
    if (!params->valid_deltamax || !params->deltamax_grid) {
        return DELTAMAX_DEFAULT;
    }
    
    int grid_idx = lat_idx * params->nlon + lon_idx;
    float deltamax_val = params->deltamax_grid[grid_idx];
    
    if (isnan(deltamax_val) || deltamax_val <= 0) {
        return DELTAMAX_DEFAULT;
    }
    
    // Opzionale: applica fattori di scala ragionevoli
    // Per CHL, i valori possono essere molto piccoli, quindi potremmo scalare
    double scaled_deltamax = (double)deltamax_val;
    
    // Limita a valori ragionevoli per evitare problemi numerici
    if (scaled_deltamax < 0.01) scaled_deltamax = 0.01;
    if (scaled_deltamax > 10.0) scaled_deltamax = 10.0;
    
    return scaled_deltamax;
}

// ===== Genera lista di punti decimati spazialmente =====
int* generate_decimated_points(int nlat, int nlon, int* total_decimated_points) {
    // Calcola quanti punti decimati avremo
    int decimated_nlat = (nlat + SPATIAL_DECIMATION - 1) / SPATIAL_DECIMATION;  // Ceil division
    int decimated_nlon = (nlon + SPATIAL_DECIMATION - 1) / SPATIAL_DECIMATION;
    *total_decimated_points = decimated_nlat * decimated_nlon;
    
    // Alloca array per gli indici dei punti decimati (formato: lat_idx * nlon + lon_idx)
    int* decimated_indices = (int*)malloc(*total_decimated_points * sizeof(int));
    
    int point_idx = 0;
    for (int lat_dec = 0; lat_dec < decimated_nlat; lat_dec++) {
        for (int lon_dec = 0; lon_dec < decimated_nlon; lon_dec++) {
            int actual_lat_idx = lat_dec * SPATIAL_DECIMATION;
            int actual_lon_idx = lon_dec * SPATIAL_DECIMATION;
            
            // Assicurati di non superare i limiti
            if (actual_lat_idx >= nlat) actual_lat_idx = nlat - 1;
            if (actual_lon_idx >= nlon) actual_lon_idx = nlon - 1;
            
            decimated_indices[point_idx] = actual_lat_idx * nlon + actual_lon_idx;
            point_idx++;
        }
    }
    
    return decimated_indices;
}

// ===== Creazione embedding da serie 1D con gestione NaN =====
VecEmb* create_embedding(const double* series, int N, int* emb_length) {
    int max_index = N - (EMB_DIM - 1) * TAU_EMB;
    if(max_index <= 0) {
        *emb_length = 0;
        return NULL;
    }
    
    VecEmb* embedded = (VecEmb*)malloc(max_index * sizeof(VecEmb));
    int valid_count = 0;
    
    for(int i = 0; i < max_index; i++) {
        int has_nan_flag = 0;
        
        // Costruisci il vettore e verifica NaN
        for(int dim = 0; dim < EMB_DIM; dim++) {
            int idx = i + dim * TAU_EMB;
            embedded[i].x[dim] = series[idx];
            if(isnan(series[idx])) {
                has_nan_flag = 1;
            }
        }
        
        embedded[i].valid = !has_nan_flag;
        if(embedded[i].valid) {
            valid_count++;
        }
    }
    
    *emb_length = max_index;
    
    // Se troppi pochi vettori validi, libera tutto
    if(valid_count < T_STEPS + 10) {
        free(embedded);
        *emb_length = 0;
        return NULL;
    }
    
    return embedded;
}

// ===== Lyapunov di Paladin con DELTAMAX adattivo =====
double lyapunov_paladin_adaptive(const double* series, int N, double deltamax_local) {
    if(N < 100) return NAN;
    
    // Crea embedding mantenendo i NaN
    int N_emb;
    VecEmb* series_emb = create_embedding(series, N, &N_emb);
    if(!series_emb || N_emb < T_STEPS + 10) {
        free(series_emb);
        return NAN;
    }
    
    double soglia_min = DELTA0 * (1 - TOL);
    double soglia_max = DELTA0 * (1 + TOL);
    
    // Trova ricorrenze considerando solo vettori validi
    int max_pairs = 50;
    int n_pairs = 0;
    int* t1_list = (int*)malloc(max_pairs * sizeof(int));
    int* t2_list = (int*)malloc(max_pairs * sizeof(int));
    int last_t1_used = -T_STEPS;
    
    for(int t1 = 0; t1 < N_emb - T_STEPS; t1++) {
        if(!series_emb[t1].valid) continue;
        if(t1 - last_t1_used < T_STEPS) continue;
        
        for(int t2 = t1 + T_STEPS; t2 < N_emb; t2++) {
            if(!series_emb[t2].valid) continue;
            
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
    
    // Calcola tempo di predicibilità con DELTAMAX adattivo
    int tau_count = 0;
    double tau_sum = 0, d0_sum = 0;
    
    for(int p = 0; p < n_pairs; p++) {
        int t1 = t1_list[p];
        int t2 = t2_list[p];
        int max_len = N_emb - (t1 > t2 ? t1 : t2);
        
        int counting = 0;
        int tau_current = 0;
        double d0_current = 0;
        int last_valid_idx1 = t1;
        int last_valid_idx2 = t2;
        double last_valid_d = 0;
        int gap_steps = 0;
        
        for(int i = 0; i < max_len; i++) {
            int idx1 = t1 + i;
            int idx2 = t2 + i;
            
            int valid1 = series_emb[idx1].valid;
            int valid2 = series_emb[idx2].valid;
            
            if(valid1 && valid2) {
                double d = norm_vec(&series_emb[idx1], &series_emb[idx2]);
                
                if(gap_steps > 0) {
                    VecEmb last_valid1 = series_emb[last_valid_idx1];
                    VecEmb last_valid2 = series_emb[last_valid_idx2];
                    double d_after_gap = norm_vec(&last_valid1, &last_valid2);
                    
                    if(d_after_gap <= deltamax_local) {  // **USA DELTAMAX ADATTIVO**
                        tau_current += gap_steps;
                        last_valid_d = d;
                        gap_steps = 0;
                    } else {
                        if(counting) {
                            tau_sum += tau_current;
                            d0_sum += d0_current;
                            tau_count++;
                            counting = 0;
                        }
                        gap_steps = 0;
                        continue;
                    }
                }
                
                if(!counting) {
                    if(d >= soglia_min && d <= soglia_max) {
                        counting = 1;
                        tau_current = 0;
                        d0_current = d;
                        last_valid_d = d;
                        last_valid_idx1 = idx1;
                        last_valid_idx2 = idx2;
                    }
                } else {
                    tau_current++;
                    if(d >= deltamax_local) {  // **USA DELTAMAX ADATTIVO**
                        tau_sum += tau_current;
                        d0_sum += d0_current;
                        tau_count++;
                        counting = 0;
                        gap_steps = 0;
                    } else {
                        last_valid_d = d;
                        last_valid_idx1 = idx1;
                        last_valid_idx2 = idx2;
                    }
                }
            } 
            else if(counting) {
                gap_steps++;
            }
        }
        
        if(counting && gap_steps == 0) {
            tau_sum += tau_current;
            d0_sum += d0_current;
            tau_count++;
        }
    }
    
    free(series_emb);
    free(t1_list);
    free(t2_list);
    
    if(tau_count == 0) return NAN;
    
    double mean_tau = tau_sum / tau_count;
    double mean_d0 = d0_sum / tau_count;
    double Tcrit = mean_tau * DT;
    double Tpred = Tcrit * log(deltamax_local / DELTA0) / log(deltamax_local / mean_d0);
    
    return Tpred;
}

// ===== Estrai serie temporale dal netCDF mantenendo i NaN =====
double* extract_series_with_nan(const float* input, int N, int* valid_count) {
    double* output = (double*)malloc(N * sizeof(double));
    *valid_count = 0;
    
    for(int i = 0; i < N; i++) {
        if(isnan(input[i]) || input[i] <= 0 || input[i] >= 100) {
            output[i] = NAN;
        } else {
            output[i] = (double)input[i];
            (*valid_count)++;
        }
    }
    
    return output;
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
    
    // Ottieni dimensioni
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
    
    // Leggi coordinate
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

// ===== Leggi un punto specifico dal netCDF =====
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

// ===== MAIN MPI =====
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc < 3 || argc > 4) {
        if (rank == 0) {
            printf("Uso: mpirun -np N %s <input_file.nc> <varname> [deltamax_file.nc]\n", argv[0]);
            printf("Esempio 1: mpirun -np 4 %s regional_dataset_CHL_10_20_-30_-15.nc CHL\n", argv[0]);
            printf("Esempio 2: mpirun -np 4 %s regional_dataset_CHL_10_20_-30_-15.nc CHL regional_dataset_CHL_10_20_-30_-15_mean.nc\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    const char* filename = argv[1];
    const char* varname = argv[2];
    const char* deltamax_filename = (argc == 4) ? argv[3] : NULL;
    
    SpatialParams params = {0};
    int total_decimated_points;
    int* decimated_indices = NULL;
    
    // Rank 0 legge le informazioni del file
    if (rank == 0) {
        printf("\n=== Lyapunov Paladin MPI with Adaptive DELTAMAX and Spatial Decimation ===\n");
        printf("Input file: %s\n", filename);
        printf("Variabile: %s\n", varname);
        if (deltamax_filename) {
            printf("DELTAMAX file: %s\n", deltamax_filename);
        } else {
            printf("DELTAMAX: costante = %f\n", DELTAMAX_DEFAULT);
        }
        printf("Decimazione spaziale: 1/%d\n", SPATIAL_DECIMATION);
        printf("Processi MPI: %d\n", size);
        
        if (!get_netcdf_info(filename, varname, &params.nlat, &params.nlon, &params.ntime, 
                            &params.lats, &params.lons)) {
            printf("ERRORE: Impossibile leggere il file netCDF\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("Dimensioni originali: lat=%d, lon=%d, time=%d\n", params.nlat, params.nlon, params.ntime);
        printf("Range lat: [%.2f, %.2f]\n", params.lats[0], params.lats[params.nlat-1]);
        printf("Range lon: [%.2f, %.2f]\n", params.lons[0], params.lons[params.nlon-1]);
        
        // Genera lista di punti decimati
        decimated_indices = generate_decimated_points(params.nlat, params.nlon, &total_decimated_points);
        printf("Punti dopo decimazione: %d (rispetto a %d originali)\n", 
               total_decimated_points, params.nlat * params.nlon);
        
        // Carica griglia DELTAMAX se fornita
        if (deltamax_filename) {
            load_deltamax_grid(deltamax_filename, &params);
        }
        
        printf("================================\n\n");
    }
    
    // Broadcast parametri a tutti i rank
    MPI_Bcast(&params.nlat, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.nlon, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.ntime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.valid_deltamax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_decimated_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Alloca e distribuisce coordinate
    if (rank != 0) {
        params.lats = (double*)malloc(params.nlat * sizeof(double));
        params.lons = (double*)malloc(params.nlon * sizeof(double));
        decimated_indices = (int*)malloc(total_decimated_points * sizeof(int));
    }
    
    MPI_Bcast(params.lats, params.nlat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(params.lons, params.nlon, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(decimated_indices, total_decimated_points, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribuisci griglia DELTAMAX
    if (params.valid_deltamax) {
        if (rank != 0) {
            params.deltamax_grid = (float*)malloc(params.nlat * params.nlon * sizeof(float));
        }
        MPI_Bcast(params.deltamax_grid, params.nlat * params.nlon, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    
    // Calcola distribuzione dei punti decimati
    int points_per_rank = total_decimated_points / size;
    int remainder = total_decimated_points % size;
    
    int start_idx, local_count;
    if (rank < remainder) {
        start_idx = rank * (points_per_rank + 1);
        local_count = points_per_rank + 1;
    } else {
        start_idx = rank * points_per_rank + remainder;
        local_count = points_per_rank;
    }
    
    printf("[Rank %d] Processando punti decimati %d - %d (totale: %d)\n", 
           rank, start_idx, start_idx + local_count - 1, local_count);
    
    // Alloca risultati locali
    double* local_paladin = (double*)malloc(local_count * sizeof(double));
    double* local_deltamax_used = (double*)malloc(local_count * sizeof(double));
    double* local_lats = (double*)malloc(local_count * sizeof(double));
    double* local_lons = (double*)malloc(local_count * sizeof(double));
    
    // Processa ogni punto decimato assegnato a questo rank
    for (int i = 0; i < local_count; i++) {
        int decimated_idx = start_idx + i;
        int global_idx = decimated_indices[decimated_idx];
        int lat_idx = global_idx / params.nlon;
        int lon_idx = global_idx % params.nlon;
        
        local_lats[i] = params.lats[lat_idx];
        local_lons[i] = params.lons[lon_idx];
        
        if (rank == 0 && i % 20 == 0) {
            printf("[Rank 0] Progresso: %d/%d (lat=%d, lon=%d)\n", i, local_count, lat_idx, lon_idx);
        }
        
        // Ottieni DELTAMAX per questo punto
        double deltamax_local = get_deltamax_for_point(&params, lat_idx, lon_idx);
        local_deltamax_used[i] = deltamax_local;
        
        // Leggi serie dal netCDF
        float* raw_data = read_netcdf_point(filename, varname, lat_idx, lon_idx, params.ntime);
        
        // Estrai serie mantenendo i NaN
        int valid_count;
        double* series = extract_series_with_nan(raw_data, params.ntime, &valid_count);
        free(raw_data);
        
        // Calcola Lyapunov con DELTAMAX adattivo
        if (valid_count >= 100) {
            local_paladin[i] = lyapunov_paladin_adaptive(series, params.ntime, deltamax_local);
        } else {
            local_paladin[i] = NAN;
        }
        
        free(series);
    }
    
    printf("[Rank %d] Calcolo completato!\n", rank);
    
    // Prepara buffer per gather
    int* recvcounts = NULL;
    int* displs = NULL;
    double* global_paladin = NULL;
    double* global_deltamax_used = NULL;
    double* global_lats = NULL;
    double* global_lons = NULL;
    
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
        
        global_paladin = (double*)malloc(total_decimated_points * sizeof(double));
        global_deltamax_used = (double*)malloc(total_decimated_points * sizeof(double));
        global_lats = (double*)malloc(total_decimated_points * sizeof(double));
        global_lons = (double*)malloc(total_decimated_points * sizeof(double));
    }
    
    // Raccogli risultati
    MPI_Gatherv(local_paladin, local_count, MPI_DOUBLE,
                global_paladin, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(local_deltamax_used, local_count, MPI_DOUBLE,
                global_deltamax_used, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(local_lats, local_count, MPI_DOUBLE,
                global_lats, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(local_lons, local_count, MPI_DOUBLE,
                global_lons, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    // Rank 0 salva i risultati
    if (rank == 0) {
        FILE* f = fopen("lyapunov_paladin_adaptive_decimated.dat", "w");
        fprintf(f, "# latitude\tlongitude\tlyapunov_paladin\tdeltamax_used\n");
        
        int valid_count = 0;
        double sum = 0, sum2 = 0;
        double min_val = 1e30, max_val = -1e30;
        double deltamax_sum = 0, deltamax_min = 1e30, deltamax_max = -1e30;
        
        for (int i = 0; i < total_decimated_points; i++) {
            fprintf(f, "%f\t%f\t%f\t%f\n", 
                    global_lats[i], global_lons[i], 
                    global_paladin[i], global_deltamax_used[i]);
            
            if (!isnan(global_paladin[i])) {
                valid_count++;
                sum += global_paladin[i];
                sum2 += global_paladin[i] * global_paladin[i];
                if (global_paladin[i] < min_val) min_val = global_paladin[i];
                if (global_paladin[i] > max_val) max_val = global_paladin[i];
                
                deltamax_sum += global_deltamax_used[i];
                if (global_deltamax_used[i] < deltamax_min) deltamax_min = global_deltamax_used[i];
                if (global_deltamax_used[i] > deltamax_max) deltamax_max = global_deltamax_used[i];
            }
        }
        fclose(f);
        
        printf("\n=== RISULTATI ===\n");
        printf("File salvato: lyapunov_paladin_adaptive_decimated.dat\n");
        printf("Punti decimati totali: %d\n", total_decimated_points);
        printf("Punti validi: %d (%.1f%%)\n", 
               valid_count, 100.0 * valid_count / total_decimated_points);
        printf("Decimazione: 1/%d (da %d punti originali)\n", 
               SPATIAL_DECIMATION, params.nlat * params.nlon);
        
        if (valid_count > 0) {
            double mean = sum / valid_count;
            double variance = (sum2 / valid_count) - (mean * mean);
            double std = sqrt(variance > 0 ? variance : 0);
            
            printf("\nLyapunov Paladin:\n");
            printf("  Media: %f\n", mean);
            printf("  Std: %f\n", std);
            printf("  Min: %f\n", min_val);
            printf("  Max: %f\n", max_val);
            
            printf("\nDELTAMAX utilizzato:\n");
            printf("  Media: %f\n", deltamax_sum / valid_count);
            printf("  Min: %f\n", deltamax_min);
            printf("  Max: %f\n", deltamax_max);
        }
        
        free(recvcounts);
        free(displs);
        free(global_paladin);
        free(global_deltamax_used);
        free(global_lats);
        free(global_lons);
    }
    
    // Libera memoria
    free(local_paladin);
    free(local_deltamax_used);
    free(local_lats);
    free(local_lons);
    free(params.lats);
    free(params.lons);
    free(decimated_indices);
    if (params.deltamax_grid) free(params.deltamax_grid);
    
    MPI_Finalize();
    return 0;
}
