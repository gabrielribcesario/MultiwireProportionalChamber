#ifndef UTILITY_HH
#define UTILITY_HH

#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <Garfield/Random.hh>

// Creates the directory that will store the results.
std::string InitializeDir(std::string base_dir){
    using namespace std;

    // Gets the current date.
    time_t t0 = time(nullptr);
    tm tm0 = *localtime(&t0);

    // Formats the date into dd/mm/yyyy format.
    stringstream ss_date;
    ss_date << tm0.tm_mday << "-" << tm0.tm_mon + 1 << "-" << tm0.tm_year + 1900;
    const string date = ss_date.str();

    // Each time the function is called a new path/date/Results/Simulation folder with incremental suffix is created.
    string tmp = base_dir + "/" + date + "/Simulation_";
    string results_dir;
    bool status = true;
    int i = 0; 
    do {
        results_dir = tmp + to_string(++i);
        status = !filesystem::create_directories(results_dir);
    }
    while(status);

    return results_dir;
}

// Samples a uniform distribution X ~ U(a,b).
double Uniform(double a, double b) { return (b - a) * Garfield::RndmUniform() + a; }

#endif