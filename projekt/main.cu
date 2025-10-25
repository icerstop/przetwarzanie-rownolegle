// Jeden plik: SEQ + OpenMP + CUDA TSP na tych samych podzbiorach

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <cuda_runtime.h>

// --------------------------------------------------
// Struktury i funkcje wsp�lne
// --------------------------------------------------
typedef struct {
    double lat, lon;
    char name[64];  // Dodane pole dla nazwy miasta
} City;

static inline double to_rad(double deg) {
    return deg * M_PI / 180.0;
}

static inline double haversine(const City *a, const City *b) {
    double dlat = to_rad(b->lat - a->lat);
    double dlon = to_rad(b->lon - a->lon);
    double rlat1 = to_rad(a->lat), rlat2 = to_rad(b->lat);
    double sdlat = sin(dlat/2), sdlon = sin(dlon/2);
    double h = sdlat*sdlat + sdlon*sdlon * cos(rlat1)*cos(rlat2);
    return 2 * 6371.0 * asin(sqrt(h));
}

void load_cities(const char *fname, City **out, int *n) {
    FILE *f = fopen(fname, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file '%s'\n", fname);
        exit(1);
    }
    char line[256];
    // Pomijamy nag��wek
    fgets(line, sizeof(line), f);
    int count = 0;
    while (fgets(line, sizeof(line), f)){
        char *tok = strtok(line, ",");
        tok = strtok(NULL, ",");
        if (!tok) break;
        tok = strtok(NULL, ",");
        if (!tok) break;
        count++;
    }
    rewind(f);
    fgets(line, sizeof(line), f);

    City *arr = (City*)malloc(count * sizeof(City));
    if (!arr){
        fprintf(stderr, "Error: Cannot allocate memory for %d cities\n", count);
        fclose(f);
        exit(1);
    }
    int i = 0;
    while (i < count && fgets(line, sizeof(line), f)) {
        char *tok = strtok(line, ",");        // nazwa miasta
        if (tok) {
            strncpy(arr[i].name, tok, sizeof(arr[i].name) - 1);
            arr[i].name[sizeof(arr[i].name) - 1] = '\0';
        } else {
            snprintf(arr[i].name, sizeof(arr[i].name), "City_%d", i);
        }

        tok = strtok(NULL, ",");
        if (!tok) break;
        arr[i].lat = atof(tok);
        tok = strtok(NULL, ",");
        if (!tok) break;
        arr[i].lon = atof(tok);
        i++;
    }
    fclose(f);
    *out = arr;
    *n = i;
}

void compute_dist(const City *cities, double *dist, int N) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            dist[i*N + j] = haversine(&cities[i], &cities[j]);
}

// --------------------------------------------------
// Funkcja do obliczania dystansu trasy
// --------------------------------------------------
double calculate_tour_distance(const int *tour, const double *dist, int N) {
    double total_distance = 0.0;
    for (int i = 0; i < N; i++) {
        int from = tour[i];
        int to = (i + 1 < N) ? tour[i + 1] : tour[0];  // Powr�t do pocz�tku
        total_distance += dist[from * N + to];
    }
    return total_distance;
}

// --------------------------------------------------
// Funkcja do wy�wietlania trasy
// --------------------------------------------------
void print_route(const int *tour, const City *cities, const double *dist, int N, const char *algorithm) {
    double distance = calculate_tour_distance(tour, dist, N);

    if (N <= 25) {  // Szczeg�owa trasa tylko dla ma�ych instancji
        printf("%s Route: ", algorithm);
        for (int i = 0; i < N; i++) {
            printf("%s", cities[tour[i]].name);
            if (i < N - 1) printf(" -> ");
        }
        printf(" -> %s (%.2f km)\n", cities[tour[0]].name, distance);
    } else {  // Tylko dystans dla wi�kszych instancji
        printf("%s Route distance: %.2f km\n", algorithm, distance);
    }
    fflush(stdout);
}

// --------------------------------------------------
// SEQ: Nearest Neighbor + 2-opt
// --------------------------------------------------
void nearest_neighbor(int *tour, const double *dist, int N) {
    int *used = (int*)calloc(N, sizeof(int));
    tour[0] = 0; used[0] = 1;
    for (int i = 1; i < N; i++) {
        int prev = tour[i-1], best = -1;
        double bestd = 1e9;
        for (int j = 0; j < N; j++) {
            if (!used[j] && dist[prev*N + j] < bestd) {
                bestd = dist[prev*N + j];
                best = j;
            }
        }
        tour[i] = best;
        used[best] = 1;
    }
    free(used);
}

void two_opt_seq(int *tour, const double *dist, int N) {
    int improved = 1;
    while (improved) {
        improved = 0;
        for (int i = 1; i < N - 1 && !improved; i++) {
            for (int j = i + 1; j < N; j++) {
                int a = tour[i-1], b = tour[i],
                    c = tour[j],   d = (j+1<N ? tour[j+1] : tour[0]);
                double delta = dist[a*N + c]
                             + dist[b*N + d]
                             - dist[a*N + b]
                             - dist[c*N + d];
                if (delta < -1e-9) {
                    for (int x = i, y = j; x < y; x++, y--) {
                        int tmp = tour[x]; tour[x] = tour[y]; tour[y] = tmp;
                    }
                    improved = 1;
                    break;
                }
            }
        }
    }
}

// --------------------------------------------------
// OpenMP: r�wnoleg�e szukanie najlepszego swapu
// --------------------------------------------------
void two_opt_omp(int *tour, const double *dist, int N) {
    int improved = 1;

    while (improved) {
        improved = 0;
        int best_i = 0, best_j = 0;
        double best_delta = 0.0;

        // r�wnoleg�a sekcja do wyszukania najlepszego lokalnego ulepszenia
        #pragma omp parallel default(none) shared(tour, dist, N, best_delta, best_i, best_j)
        {
            int loc_i = 0, loc_j = 0;
            double loc_best = 0.0;

            // r�wnoleg�e roz�o�enie iteracji z pomini�ciem synchronizacji ko�ca p�tli
            #pragma omp for nowait
            for (int i = 1; i < N - 1; i++) {
                for (int j = i + 1; j < N; j++) {
                    int a = tour[i - 1], b = tour[i];
                    int c = tour[j], d = (j + 1 < N ? tour[j + 1] : tour[0]);

                    double delta = dist[a * N + c] + dist[b * N + d]
                                 - dist[a * N + b] - dist[c * N + d];

                    if (delta < loc_best) {
                        loc_best = delta;
                        loc_i = i;
                        loc_j = j;
                    }
                }
            }

            // sekcja krytyczna � aktualizacja najlepszej znanej poprawy
            #pragma omp critical
            {
                if (loc_best < best_delta) {
                    best_delta = loc_best;
                    best_i = loc_i;
                    best_j = loc_j;
                }
            }
        }

        // Zastosowanie najlepszego znalezionego ulepszenia (je�li wyst�puje)
        if (best_delta < -1e-9) {
            for (int x = best_i, y = best_j; x < y; x++, y--) {
                int tmp = tour[x];
                tour[x] = tour[y];
                tour[y] = tmp;
            }
            improved = 1;
        }
    }
}
// --------------------------------------------------
// CUDA: atomowa operacja min dla double
// --------------------------------------------------
__device__ double atomicMinDouble(double *addr, double val) {
    unsigned long long *ptr = (unsigned long long*)addr;
    unsigned long long old = *ptr, assumed;
    do {
        assumed = old;
        if (__longlong_as_double(assumed) <= val) break;
        old = atomicCAS(ptr, assumed, __double_as_longlong(val));
    } while (assumed != old);
    return __longlong_as_double(old);
}

// --------------------------------------------------
// CUDA kernel: wyszukaj najlepszy swap
// --------------------------------------------------
__global__ void find_best_swap(int *tour, double *dist, int N,
                               int *best_i, int *best_j, double *best_delta) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = (N-1)*N/2;
    if (idx >= total) return;

    // zamapuj idx -> (i,j)
    int i = 1, rem = idx;
    while (rem >= N - i) { rem -= (N - i); i++; }
    int j = i + 1 + rem;

    int a = tour[i-1], b = tour[i],
        c = tour[j],   d = (j+1<N ? tour[j+1] : tour[0]);
    double delta = dist[a*N + c]
                 + dist[b*N + d]
                 - dist[a*N + b]
                 - dist[c*N + d];

    // atomowo zaktualizuj najlepsz� popraw�
    double prev = atomicMinDouble(best_delta, delta);
    if (prev > delta) {
        atomicExch(best_i, i);
        atomicExch(best_j, j);
    }
}

// --------------------------------------------------
// GPU: 2-opt z kernelem
// --------------------------------------------------
void two_opt_cuda(int *tour, double *dist, int N) {
    int *d_tour, *d_best_i, *d_best_j;
    double *d_dist, *d_best_delta;
    int total_pairs = (N-1)*N/2;

    cudaMalloc(&d_tour,       N * sizeof(int));
    cudaMalloc(&d_dist,       N * N * sizeof(double));
    cudaMalloc(&d_best_i,     sizeof(int));
    cudaMalloc(&d_best_j,     sizeof(int));
    cudaMalloc(&d_best_delta, sizeof(double));

    cudaMemcpy(d_tour, tour, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dist, dist, N * N * sizeof(double), cudaMemcpyHostToDevice);

    int improved = 1;
    while (improved) {
        improved = 0;
        double h_delta = 0.0;
        int h_i = 0, h_j = 0;
        cudaMemcpy(d_best_delta, &h_delta, sizeof(double), cudaMemcpyHostToDevice);

        int threads = 1024;
        int blocks  = (total_pairs + threads - 1) / threads;
        find_best_swap<<<blocks, threads>>>(d_tour, d_dist, N,
                                            d_best_i, d_best_j, d_best_delta);
        cudaDeviceSynchronize();

        cudaMemcpy(&h_delta,   d_best_delta, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&h_i,       d_best_i,     sizeof(int),    cudaMemcpyDeviceToHost);
        cudaMemcpy(&h_j,       d_best_j,     sizeof(int),    cudaMemcpyDeviceToHost);

        if (h_delta < -1e-9) {
            // zastosuj swap na CPU
            for (int x = h_i, y = h_j; x < y; x++, y--) {
                int tmp = tour[x]; tour[x] = tour[y]; tour[y] = tmp;
            }
            // zsynchronizuj z GPU i powt�rz
            cudaMemcpy(d_tour, tour, N * sizeof(int), cudaMemcpyHostToDevice);
            improved = 1;
        }
    }

    cudaMemcpy(tour, d_tour, N * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_tour);       cudaFree(d_dist);
    cudaFree(d_best_i);     cudaFree(d_best_j);
    cudaFree(d_best_delta);
}

// --------------------------------------------------
// main(): argumenty z linii komend
// --------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <cities_file> <size1> <size2> ... <sizeN>\n", argv[0]);
        fprintf(stderr, "Example: %s cities.csv 10 25 50 150 250\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    int num_sizes = argc - 2;
    int *sizes = (int*)malloc(num_sizes * sizeof(int));

    // Parsowanie rozmiar�w z argument�w
    for (int i = 0; i < num_sizes; i++) {
        sizes[i] = atoi(argv[i + 2]);
        if (sizes[i] <= 0) {
            fprintf(stderr, "Error: Size must be positive integer, got '%s'\n", argv[i + 2]);
            free(sizes);
            return 1;
        }
    }

    City *cities;
    int fullN;
    load_cities(filename, &cities, &fullN);
    printf("Loaded %d cities from '%s'\n", fullN, filename);

    // Sprawdzenie czy wszystkie rozmiary s� mo�liwe
    for (int i = 0; i < num_sizes; i++) {
        if (sizes[i] > fullN) {
            fprintf(stderr, "Error: Requested size %d is larger than available cities (%d)\n",
                    sizes[i], fullN);
            free(sizes);
            free(cities);
            return 1;
        }
    }

    cudaFree(0);
    int max_procs = omp_get_num_procs();
    omp_set_num_threads(max_procs);

    #pragma omp parallel
    {}

    {
        int dummyTour[2] = {0,1};
        double dummyDist[4] = {0.0,0.0,0.0,0.0};
        two_opt_omp(dummyTour, dummyDist, 2);
    }

    const int runs = 5;
    FILE *f;

    // nag��wek
    f = fopen("results.csv", "w");
    fprintf(f, "alg,size,time\n");
    fclose(f);

    srand(time(NULL));
    printf("Starting benchmark with %d runs for each size and algorithm\n", runs);

    for (int r = 0; r < runs; r++) {
        printf("Run %d/%d\n", r + 1, runs);
        for (int si = 0; si < num_sizes; si++) {
            int subN = sizes[si];

            // wsp�lny podzbi�r
            int *idx    = (int*)malloc(subN * sizeof(int));
            int *picked = (int*)calloc(fullN, sizeof(int));
            for (int i = 0; i < subN; i++) {
                int x;
                do { x = rand() % fullN; } while (picked[x]);
                picked[x] = 1;
                idx[i]    = x;
            }
            free(picked);

            // budujemy macierz odleg�o�ci
            City   *sub  = (City*)  malloc(subN * sizeof(City));
            double *dist = (double*)malloc(subN * subN * sizeof(double));
            for (int i = 0; i < subN; i++) sub[i] = cities[idx[i]];
            compute_dist(sub, dist, subN);

            int *tour = (int*)malloc(subN * sizeof(int));
            struct timespec t1, t2;
            double dt;

            // SEQ
            clock_gettime(CLOCK_MONOTONIC, &t1);
            nearest_neighbor(tour, dist, subN);
            two_opt_seq(tour, dist, subN);
            clock_gettime(CLOCK_MONOTONIC, &t2);
            dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
            f = fopen("results.csv","a");
            fprintf(f, "SEQ,%d,%.6f\n", subN, dt);
            fclose(f);
            printf("Done SEQ: %d cities, %.6f s\n", subN, dt);
            print_route(tour, sub, dist, subN, "SEQ");  // Wy�wietl tras� dla ma�ych instancji
            fflush(stdout);

            // OpenMP
            clock_gettime(CLOCK_MONOTONIC, &t1);
            nearest_neighbor(tour, dist, subN);
            two_opt_omp(tour, dist, subN);
            clock_gettime(CLOCK_MONOTONIC, &t2);
            dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
            f = fopen("results.csv","a");
            fprintf(f, "OMP,%d,%.6f\n", subN, dt);
            fclose(f);
            printf("Done OMP: %d cities, %.6f s\n", subN, dt);
            print_route(tour, sub, dist, subN, "OMP");  // Wy�wietl tras� dla ma�ych instancji
            fflush(stdout);

            // CUDA
            clock_gettime(CLOCK_MONOTONIC, &t1);
            nearest_neighbor(tour, dist, subN);
            two_opt_cuda(tour, dist, subN);
            clock_gettime(CLOCK_MONOTONIC, &t2);
            dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
            f = fopen("results.csv","a");
            fprintf(f, "CUDA,%d,%.6f\n", subN, dt);
            fclose(f);
            printf("Done CUDA: %d cities, %.6f s\n", subN, dt);
            print_route(tour, sub, dist, subN, "CUDA");  // Wy�wietl tras� dla ma�ych instancji
            fflush(stdout);

            free(idx); free(sub); free(dist); free(tour);
        }
    }

    free(sizes);
    free(cities);
    return 0;
