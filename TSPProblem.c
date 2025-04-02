#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <time.h>
#include <float.h>

#define PI 3.14159265358979323846
#define EARTH_RADIUS_KM 6371.0

typedef struct {
    char name[50];
    double lat, lon;
    int priority;
} City;

// Pomocnicze funkcje do sortowania (potrzebne przy bruteForceTSP)
int compare_cities(const void* a, const void* b) {
    const City* city_a = (const City*)a;
    const City* city_b = (const City*)b;
    return strcmp(city_a->name, city_b->name);
}

bool cities_equal(City city1, City city2) {
    return strcmp(city1.name, city2.name) == 0;
}

// Funkcja generująca losowe miasta
City* generateRandomCities(int n) {
    City* cities = (City*)malloc(n * sizeof(City));
    
    // Używamy aktualnego czasu jako ziarna dla generatora
    srand((unsigned int)time(NULL));
    
    for (int i = 0; i < n; ++i) {
        char buffer[50];
        sprintf(buffer, "City%d", i);
        strcpy(cities[i].name, buffer);
        
        // Zakres szerokości geograficznej Polski: 49.0 - 55.0
        cities[i].lat = 49.0 + ((double)rand() / RAND_MAX) * 6.0;
        
        // Zakres długości geograficznej Polski: 14.0 - 24.0
        cities[i].lon = 14.0 + ((double)rand() / RAND_MAX) * 10.0;
        
        // Priorytet: 1, 2 lub 3
        cities[i].priority = 1 + rand() % 3;
    }
    
    return cities;
}

double toRadians(double degree) {
    return degree * (PI / 180.0);
}

double calculateDistance(double lat1, double lon1, double lat2, double lon2) {
    double lat1Rad = toRadians(lat1);
    double lat2Rad = toRadians(lat2);
    double deltaLat = toRadians(lat2 - lat1);
    double deltaLon = toRadians(lon2 - lon1);
    
    double a = sin(deltaLat / 2) * sin(deltaLat / 2) +
               cos(lat1Rad) * cos(lat2Rad) *
               sin(deltaLon / 2) * sin(deltaLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    return EARTH_RADIUS_KM * c;
}

int findStartingCityIndex(const City* cities, int size) {
    int highestPriority = INT_MAX;
    int* highestPriorityCities = (int*)malloc(size * sizeof(int));
    int highestPriorityCitiesCount = 0;
    
    // Znajdź miasta o najwyższym priorytecie
    for (int i = 0; i < size; ++i) {
        if (cities[i].priority < highestPriority) {
            highestPriority = cities[i].priority;
            highestPriorityCities[0] = i;
            highestPriorityCitiesCount = 1;
        }
        else if (cities[i].priority == highestPriority) {
            highestPriorityCities[highestPriorityCitiesCount++] = i;
        }
    }
    
    // Znajdź miasto o najwyższym priorytecie, które jest najbliżej centralnego punktu geograficznego Polski
    double centralLat = 52.215933;
    double centralLon = 19.134422;
    double minDistance = DBL_MAX;
    int startIndex = -1;
    
    for (int i = 0; i < highestPriorityCitiesCount; ++i) {
        int cityIndex = highestPriorityCities[i];
        double distance = calculateDistance(centralLat, centralLon, cities[cityIndex].lat, cities[cityIndex].lon);
        if (distance < minDistance) {
            minDistance = distance;
            startIndex = cityIndex;
        }
    }
    
    free(highestPriorityCities);
    return startIndex;
}

City* greedyTSP(const City* cities, int size, int startIndex, double* totalDistance) {
    City* tour = (City*)malloc((size + 1) * sizeof(City)); // +1 dla powrotu do punktu początkowego
    bool* visited = (bool*)calloc(size, sizeof(bool)); // Inicjalizacja wszystkich na false
    
    int currentCityIndex = startIndex;
    visited[currentCityIndex] = true;
    tour[0] = cities[currentCityIndex];
    
    *totalDistance = 0.0;
    int tourSize = 1;
    
    while (tourSize < size) {
        double minDistance = DBL_MAX;
        int highestPriority = INT_MAX;
        int nextCityIndex = -1;
        
        for (int i = 0; i < size; ++i) {
            if (!visited[i]) {
                double distance = calculateDistance(
                    cities[currentCityIndex].lat, cities[currentCityIndex].lon,
                    cities[i].lat, cities[i].lon
                );
                
                if (cities[i].priority < highestPriority ||
                   (cities[i].priority == highestPriority && distance < minDistance)) {
                    highestPriority = cities[i].priority;
                    minDistance = distance;
                    nextCityIndex = i;
                }
            }
        }
        
        if (nextCityIndex == -1) break; // W przypadku gdy nie ma więcej miast do odwiedzenia
        
        currentCityIndex = nextCityIndex;
        visited[currentCityIndex] = true;
        tour[tourSize++] = cities[currentCityIndex];
        *totalDistance += minDistance; // Dodaj odległość do łącznej sumy
    }
    
    // Dodanie odległości powrotu do punktu startowego
    if (tourSize > 0) {
        *totalDistance += calculateDistance(
            tour[tourSize-1].lat, tour[tourSize-1].lon,
            cities[startIndex].lat, cities[startIndex].lon
        );
        tour[tourSize++] = cities[startIndex]; // Dodanie punktu startowego na końcu trasy
    }
    
    // Wyświetlanie trasy i łącznego dystansu
    printf("\nGreedy TSP route: ");
    for (int i = 0; i < tourSize; ++i) {
        printf("%s -> ", tour[i].name);
    }
    printf("\b\b\b   \n"); // Usuwa ostatnią strzałkę
    printf("Total distance: %.2f km\n", *totalDistance);
    
    return tour;
}

City* randomGreedyTSP(const City* cities, int size, int startIndex, double* totalDistance) {
    City* tour = (City*)malloc((size + 1) * sizeof(City)); // +1 dla powrotu do punktu początkowego
    bool* visited = (bool*)calloc(size, sizeof(bool)); // Inicjalizacja wszystkich na false
    
    int currentCityIndex = startIndex;
    visited[currentCityIndex] = true;
    tour[0] = cities[currentCityIndex];
    
    *totalDistance = 0.0;
    int tourSize = 1;
    
    while (tourSize < size) {
        int* candidates = (int*)malloc(size * sizeof(int));
        int candidatesCount = 0;
        int currentPriority = INT_MAX;
        
        // Znajdź miasta o najwyższym dostępnym priorytecie
        for (int i = 0; i < size; ++i) {
            if (!visited[i]) {
                if (cities[i].priority < currentPriority) {
                    candidatesCount = 0;
                    currentPriority = cities[i].priority;
                    candidates[candidatesCount++] = i;
                }
                else if (cities[i].priority == currentPriority) {
                    candidates[candidatesCount++] = i;
                }
            }
        }
        
        if (candidatesCount == 0) break; // Jeśli nie ma więcej miast do odwiedzenia
        
        // Wybierz losowo miasto z kandydatów
        int randomIndex = rand() % candidatesCount;
        int nextCityIndex = candidates[randomIndex];
        
        // Dodaj odległość do łącznej sumy
        *totalDistance += calculateDistance(
            cities[currentCityIndex].lat, cities[currentCityIndex].lon,
            cities[nextCityIndex].lat, cities[nextCityIndex].lon
        );
        
        currentCityIndex = nextCityIndex;
        visited[currentCityIndex] = true;
        tour[tourSize++] = cities[currentCityIndex];
        
        free(candidates);
    }
    
    // Dodanie odległości powrotu do punktu startowego
    if (tourSize > 0) {
        *totalDistance += calculateDistance(
            tour[tourSize-1].lat, tour[tourSize-1].lon,
            cities[startIndex].lat, cities[startIndex].lon
        );
        tour[tourSize++] = cities[startIndex]; // Dodanie punktu startowego na końcu trasy
    }
    
    // Wyświetlanie trasy i łącznego dystansu
    printf("\nRandom Greedy TSP route: ");
    for (int i = 0; i < tourSize; ++i) {
        printf("%s -> ", tour[i].name);
    }
    printf("\b\b\b   \n"); // Usuwa ostatnią strzałkę
    printf("Total distance: %.2f km\n\n", *totalDistance);
    
    return tour;
}

// Struktura dla miast pogrupowanych według priorytetu
typedef struct {
    int priority;
    City* cities;
    int count;
} PriorityGroup;

void swap(City* a, City* b) {
    City temp = *a;
    *a = *b;
    *b = temp;
}

bool next_permutation(City* arr, int size) {
    if (size <= 1) return false;
    
    int i = size - 2;
    while (i >= 0 && compare_cities(&arr[i], &arr[i+1]) >= 0) i--;
    
    if (i < 0) return false;
    
    int j = size - 1;
    while (compare_cities(&arr[j], &arr[i]) <= 0) j--;
    
    swap(&arr[i], &arr[j]);
    
    // Odwrócenie części tablicy od i+1 do końca
    for (int k = i+1, l = size-1; k < l; k++, l--) {
        swap(&arr[k], &arr[l]);
    }
    
    return true;
}

double calculateTotalDistance(PriorityGroup* groups, int groupCount) {
    double totalDistance = 0.0;
    City lastCityInPreviousPriority;
    bool isFirstGroup = true;
    
    for (int i = 0; i < groupCount; ++i) {
        if (groups[i].count > 0) {
            // Oblicz odległość między grupami priorytetów
            if (!isFirstGroup) {
                totalDistance += calculateDistance(
                    lastCityInPreviousPriority.lat, lastCityInPreviousPriority.lon,
                    groups[i].cities[0].lat, groups[i].cities[0].lon
                );
            }
            
            // Oblicz odległości wewnątrz grupy priorytetów
            for (int j = 0; j < groups[i].count - 1; ++j) {
                totalDistance += calculateDistance(
                    groups[i].cities[j].lat, groups[i].cities[j].lon,
                    groups[i].cities[j+1].lat, groups[i].cities[j+1].lon
                );
            }
            
            lastCityInPreviousPriority = groups[i].cities[groups[i].count - 1];
            isFirstGroup = false;
        }
    }
    
    // Dodajemy dystans powrotny do pierwszego miasta pierwszej grupy
    if (groupCount > 0 && groups[0].count > 0) {
        totalDistance += calculateDistance(
            lastCityInPreviousPriority.lat, lastCityInPreviousPriority.lon,
            groups[0].cities[0].lat, groups[0].cities[0].lon
        );
    }
    
    return totalDistance;
}

void bruteForceTSP(City* cities, int size) {
    // Grupowanie miast według priorytetu
    int maxPriority = 0;
    for (int i = 0; i < size; ++i) {
        if (cities[i].priority > maxPriority) {
            maxPriority = cities[i].priority;
        }
    }
    
    PriorityGroup* groups = (PriorityGroup*)malloc((maxPriority + 1) * sizeof(PriorityGroup));
    for (int i = 1; i <= maxPriority; ++i) {
        groups[i-1].priority = i;
        groups[i-1].count = 0;
        
        // Liczymy miasta z danym priorytetem
        for (int j = 0; j < size; ++j) {
            if (cities[j].priority == i) {
                groups[i-1].count++;
            }
        }
        
        groups[i-1].cities = (City*)malloc(groups[i-1].count * sizeof(City));
        
        // Wypełniamy tablicę miast dla danego priorytetu
        int index = 0;
        for (int j = 0; j < size; ++j) {
            if (cities[j].priority == i) {
                groups[i-1].cities[index++] = cities[j];
            }
        }
    }
    
    double minDistance = DBL_MAX;
    City* bestRoute = (City*)malloc(size * sizeof(City));
    
    // Proste przeszukiwanie wszystkich permutacji
    bool keepGoing = true;
    while (keepGoing) {
        double currentDistance = calculateTotalDistance(groups, maxPriority);
        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            
            // Zapisz najlepszą trasę
            int index = 0;
            for (int i = 0; i < maxPriority; ++i) {
                for (int j = 0; j < groups[i].count; ++j) {
                    bestRoute[index++] = groups[i].cities[j];
                }
            }
        }
        
        // Znajdź następną permutację
        keepGoing = false;
        for (int i = maxPriority - 1; i >= 0; --i) {
            if (next_permutation(groups[i].cities, groups[i].count)) {
                keepGoing = true;
                break;
            }
        }
    }
    
    printf("Brute-force TSP route: ");
    for (int i = 0; i < size; ++i) {
        printf("%s -> ", bestRoute[i].name);
    }
    printf("%s\n", bestRoute[0].name); // Powrót do punktu startowego
    printf("Total distance: %.2f km\n", minDistance);
    
    // Zwolnienie pamięci
    for (int i = 0; i < maxPriority; ++i) {
        free(groups[i].cities);
    }
    free(groups);
    free(bestRoute);
}

int main() {
    City test_cities[] = {
        {"Warsaw", 52.23, 21.01, 1},
        {"Krakow", 50.06, 19.95, 2},
        {"Gdansk", 54.35, 18.64, 2},
        {"Lodz", 51.76, 19.46, 2},
        {"Wroclaw", 51.11, 17.04, 2},
        {"Poznan", 52.4, 16.93, 2},
        {"Szczecin", 53.43, 14.56, 3},
        {"Lublin", 51.25, 22.57, 3},
        {"Katowice", 50.27, 19.04, 3},
        {"Bydgoszcz", 53.12, 18.01, 3},
        {"Bialystok", 53.13, 23.17, 3},
        {"Radom", 51.40, 21.15, 3}
    };
    
    int test_cities_size = sizeof(test_cities) / sizeof(test_cities[0]);
    
    for (int i = 0; i < test_cities_size; ++i) {
        printf("City: %s, Lat: %.2f, Lon: %.2f, Priority: %d\n", 
               test_cities[i].name, test_cities[i].lat, test_cities[i].lon, test_cities[i].priority);
    }
    
    printf("\n");
    
    // Generator losowych miast
    int n;
    printf("Enter the number of cities: ");
    scanf("%d", &n);
    City* cities = generateRandomCities(n);
    
    // Wyświetlenie wygenerowanych miast
    for (int i = 0; i < n; ++i) {
        printf("City: %s, Lat: %.2f, Lon: %.2f, Priority: %d\n", 
               cities[i].name, cities[i].lat, cities[i].lon, cities[i].priority);
    }
    
    // Algorytm zachłanny
    double totalDistanceGreedy;
    int startIndex = findStartingCityIndex(test_cities, test_cities_size);
    City* tourGreedy = greedyTSP(test_cities, test_cities_size, startIndex, &totalDistanceGreedy);
    
    // Algorytm zachłanny losowy
    double totalDistanceRandomGreedy;
    City* tourRandomGreedy = randomGreedyTSP(test_cities, test_cities_size, startIndex, &totalDistanceRandomGreedy);
    
    // Algorytm brute force
    bruteForceTSP(test_cities, test_cities_size);
    
    // Zwolnienie pamięci
    free(tourGreedy);
    free(tourRandomGreedy);
    free(cities);
    
    return 0;
}