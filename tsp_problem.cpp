#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <limits>
#include <map>
#include <chrono>
#include <algorithm>

struct City {
    std::string name;
    double lat, lon;
    int priority;
    bool operator<(const City& other) const {
        return name < other.name;
    }
    bool operator==(const City& other) const {
        return name == other.name;
    }
};

std::vector<City> generateRandomCities(int n) {
    std::vector<City> cities;

    // Użycie aktualnego czasu jako ziarna dla generatora
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    std::uniform_real_distribution<double> lat_distribution(49.0, 55.0); // Zakres szerokości geograficznej Polski
    std::uniform_real_distribution<double> lon_distribution(14.0, 24.0); // Zakres długości geograficznej Polski
    std::uniform_int_distribution<int> priority_distribution(1, 3);

    for (int i = 0; i < n; ++i) {
        City city;
        city.name = "City" + std::to_string(i);
        city.lat = lat_distribution(generator);
        city.lon = lon_distribution(generator);
        city.priority = priority_distribution(generator);

        cities.push_back(city);
    }

    return cities;
}

const double PI = 3.14159265358979323846;
const double EARTH_RADIUS_KM = 6371.0;

double toRadians(double degree) {
    return degree * (PI / 180.0);
}

double calculateDistance(double lat1, double lon1, double lat2, double lon2) {
    double lat1Rad = toRadians(lat1);
    double lat2Rad = toRadians(lat2);
    double deltaLat = toRadians(lat2 - lat1);
    double deltaLon = toRadians(lon2 - lon1);

    double a = std::sin(deltaLat / 2) * std::sin(deltaLat / 2) +
        std::cos(lat1Rad) * std::cos(lat2Rad) *
        std::sin(deltaLon / 2) * std::sin(deltaLon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    return EARTH_RADIUS_KM * c;
}

int findStartingCityIndex(const std::vector<City>& cities) {
    int highestPriority = std::numeric_limits<int>::max();
    std::vector<int> highestPriorityCities;

    // Znajdź miasta o najwyższym priorytecie
    for (int i = 0; i < cities.size(); ++i) {
        if (cities[i].priority < highestPriority) {
            highestPriority = cities[i].priority;
            highestPriorityCities.clear();
            highestPriorityCities.push_back(i);
        }
        else if (cities[i].priority == highestPriority) {
            highestPriorityCities.push_back(i);
        }
    }

    // Teraz znajdź miasto o najwyższym priorytecie, które jest najbliżej centralnego punktu geograficznego Polski
    double centralLat = 52.215933;
    double centralLon = 19.134422;
    double minDistance = std::numeric_limits<double>::max();
    int startIndex = -1;

    for (int i : highestPriorityCities) {
        double distance = calculateDistance(centralLat, centralLon, cities[i].lat, cities[i].lon);
        if (distance < minDistance) {
            minDistance = distance;
            startIndex = i;
        }
    }

    return startIndex;
}

std::vector<City> greedyTSP(const std::vector<City>& cities, int startIndex, double& totalDistance) {
    std::vector<City> tour;
    std::vector<bool> visited(cities.size(), false);
    int currentCityIndex = startIndex;
    visited[currentCityIndex] = true;
    tour.push_back(cities[currentCityIndex]);

    totalDistance = 0.0;

    while (tour.size() < cities.size()) {
        double minDistance = std::numeric_limits<double>::max();
        int highestPriority = std::numeric_limits<int>::max();
        int nextCityIndex = -1;

        for (int i = 0; i < cities.size(); ++i) {
            if (!visited[i]) {
                double distance = calculateDistance(cities[currentCityIndex].lat, cities[currentCityIndex].lon, cities[i].lat, cities[i].lon);
                if (cities[i].priority < highestPriority ||
                    (cities[i].priority == highestPriority && distance < minDistance)) {
                    highestPriority = cities[i].priority;
                    minDistance = distance;
                    nextCityIndex = i;
                }
            }
        }

        if (nextCityIndex == -1) break; // W przypadku, gdy nie ma więcej miast do odwiedzenia
        currentCityIndex = nextCityIndex;
        visited[currentCityIndex] = true;
        tour.push_back(cities[currentCityIndex]);
        totalDistance += minDistance; // Dodaj odległość do łącznej sumy
    }

    // Dodanie odległości powrotu do punktu startowego
    if (!tour.empty()) {
        totalDistance += calculateDistance(tour.back().lat, tour.back().lon, cities[startIndex].lat, cities[startIndex].lon);
        tour.push_back(cities[startIndex]); // Dodanie punktu startowego na końcu trasy
    }

    // Wyświetlanie trasy i łącznego dystansu
    std::cout << "\nGreedy TSP route: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i].name;
        if (i < tour.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
    std::cout << "Total distance: " << totalDistance << " km\n";

    return tour;
}

std::vector<City> randomGreedyTSP(const std::vector<City>& cities, int startIndex, double& totalDistance) {
    std::vector<City> tour;
    std::vector<bool> visited(cities.size(), false);
    int currentCityIndex = startIndex;
    visited[currentCityIndex] = true;
    tour.push_back(cities[currentCityIndex]);

    // Inicjowanie generatora liczb losowych z ziarnem
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    totalDistance = 0.0;

    while (tour.size() < cities.size()) {
        std::vector<int> candidates; // Kandydaci do następnego miasta
        int currentPriority = std::numeric_limits<int>::max();

        // Znajdź miasta o najwyższym dostępnym priorytecie
        for (int i = 0; i < cities.size(); ++i) {
            if (!visited[i]) {
                if (cities[i].priority < currentPriority) {
                    candidates.clear();
                    currentPriority = cities[i].priority;
                    candidates.push_back(i);
                }
                else if (cities[i].priority == currentPriority) {
                    candidates.push_back(i);
                }
            }
        }

        if (candidates.empty()) break; // Jeśli nie ma więcej miast do odwiedzenia

        // Wybierz losowo miasto z kandydatów
        std::uniform_int_distribution<int> distribution(0, candidates.size() - 1);
        int nextCityIndex = candidates[distribution(generator)];

        // Dodaj odległość do łącznej sumy
        totalDistance += calculateDistance(cities[currentCityIndex].lat, cities[currentCityIndex].lon, cities[nextCityIndex].lat, cities[nextCityIndex].lon);

        currentCityIndex = nextCityIndex;
        visited[currentCityIndex] = true;
        tour.push_back(cities[currentCityIndex]);
    }

    // Dodanie odległości powrotu do punktu startowego
    if (!tour.empty()) {
        totalDistance += calculateDistance(tour.back().lat, tour.back().lon, cities[startIndex].lat, cities[startIndex].lon);
        tour.push_back(cities[startIndex]); // Dodanie punktu startowego na końcu trasy
    }

    // Wyświetlanie trasy i łącznego dystansu
    std::cout << "\nRandom Greedy TSP route: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i].name;
        if (i < tour.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
    std::cout << "Total distance: " << totalDistance << " km\n\n";

    return tour;
}

// Poprawiona implementacja bruteForceTSP dla małych zestawów danych
void bruteForceTSP(const std::vector<City>& cities) {
    if (cities.size() > 11) {
        std::cout << "Warning: Brute force approach not recommended for more than 10 cities." << std::endl;
        std::cout << "The computation might take extremely long time. Continue anyway? (y/n): ";
        char answer;
        std::cin >> answer;
        if (answer != 'y' && answer != 'Y') {
            std::cout << "Brute force TSP cancelled." << std::endl;
            return;
        }
    }

    // Kopiujemy miasta, aby ich nie modyfikować
    std::vector<City> allCities = cities;
    
    // Upewniamy się, że nie mamy zduplikowanego miasta startowego na końcu
    if (allCities.size() > 1 && allCities.front().name == allCities.back().name) {
        allCities.pop_back();
    }
    
    // Ustalamy miasto startowe
    City startCity = allCities.front();
    
    // Usuwamy miasto startowe z listy do permutacji
    std::vector<City> citiesToPermute(allCities.begin() + 1, allCities.end());
    
    double minDistance = std::numeric_limits<double>::max();
    std::vector<City> bestRoute;
    
    // Liczymy ile permutacji będziemy sprawdzać
    unsigned long long totalPermutations = 1;
    for (int i = 1; i <= citiesToPermute.size(); i++) {
        totalPermutations *= i;
    }
    
    std::cout << "Checking " << totalPermutations << " possible routes..." << std::endl;
    
    unsigned long long currentPermutation = 0;
    
    // Generujemy wszystkie możliwe permutacje miast (oprócz pierwszego)
    do {
        std::vector<City> currentRoute;
        currentRoute.push_back(startCity); // Dodajemy miasto startowe
        
        // Dodajemy pozostałe miasta w aktualnej permutacji
        currentRoute.insert(currentRoute.end(), citiesToPermute.begin(), citiesToPermute.end());
        
        // Dodajemy miasto startowe na końcu, aby zamknąć cykl
        currentRoute.push_back(startCity);
        
        // Obliczamy długość trasy
        double currentDistance = 0.0;
        for (size_t i = 0; i < currentRoute.size() - 1; ++i) {
            currentDistance += calculateDistance(
                currentRoute[i].lat, currentRoute[i].lon,
                currentRoute[i + 1].lat, currentRoute[i + 1].lon
            );
        }
        
        // Wyświetlamy postęp co 1000 permutacji
        currentPermutation++;
        if (currentPermutation % 1000 == 0) {
            std::cout << "Checked " << currentPermutation << " of " << totalPermutations << " routes (" 
                      << (currentPermutation * 100.0 / totalPermutations) << "%)\r" << std::flush;
        }
        
        // Aktualizujemy najlepszą trasę, jeśli znaleźliśmy lepszą
        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            bestRoute = currentRoute;
        }
        
    } while (std::next_permutation(citiesToPermute.begin(), citiesToPermute.end(), 
                                   [](const City& a, const City& b) { return a.name < b.name; }));
    
    std::cout << std::endl;
    
    // Wyświetlamy najlepszą trasę
    std::cout << "Brute-force TSP route: ";
    for (size_t i = 0; i < bestRoute.size(); ++i) {
        std::cout << bestRoute[i].name;
        if (i < bestRoute.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
    std::cout << "Total distance: " << minDistance << " km" << std::endl;
}

int main() {
    std::vector<City> test_cities = {
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

    std::cout << "Pre-defined test cities:" << std::endl;
    for (const auto& city : test_cities) {
        std::cout << "City: " << city.name
            << ", Lat: " << city.lat
            << ", Lon: " << city.lon
            << ", Priority: " << city.priority << std::endl;
    }

    std::cout << std::endl;
    
    // Poproś użytkownika o wybór zestawu danych
    char choice;
    std::cout << "Do you want to use pre-defined cities or generate random ones? (p/r): ";
    std::cin >> choice;
    
    std::vector<City> cities;
    if (choice == 'p' || choice == 'P') {
        cities = test_cities;
    } else {
        // Generator losowych miast
        int n; // Liczba miast
        std::cout << "Enter the number of cities: ";
        std::cin >> n;
        cities = generateRandomCities(n);

        // Wyświetlenie wygenerowanych miast
        std::cout << "Generated random cities:" << std::endl;
        for (const auto& city : cities) {
            std::cout << "City: " << city.name
                << ", Lat: " << city.lat
                << ", Lon: " << city.lon
                << ", Priority: " << city.priority << std::endl;
        }
        std::cout << std::endl;
    }

    // Algorytm zachłanny
    double totalDistanceGreedy;
    int startIndex = findStartingCityIndex(cities);
    std::cout << "Starting city: " << cities[startIndex].name << std::endl;
    
    std::vector<City> tourGreedy = greedyTSP(cities, startIndex, totalDistanceGreedy);

    // Algorytm zachłanny losowy
    double totalDistanceRandomGreedy;
    std::vector<City> tourRandomGreedy = randomGreedyTSP(cities, startIndex, totalDistanceRandomGreedy);

    // Brute Force TSP (tylko dla małych zestawów danych)
    bruteForceTSP(cities);

    return 0;
}