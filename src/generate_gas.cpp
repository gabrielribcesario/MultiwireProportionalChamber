#include <cstdlib>
#include <filesystem>
#include <string>

#include <Garfield/MediumMagboltz.hh>
#include <Garfield/FundamentalConstants.hh>

/*
    Creates a table of electron transport parameters for a given gas mixture and range of electric fields.
*/

using namespace Garfield;

int main(int argc, char * argv[]) {
    // Number of collisions (in multiples of 10^7) over which the electron is traced by Magboltz.
    // Recommended value is nCollisions >= 10.
    const int nCollisions = 10;
    // Number of electric field points to be evaluated.
    const size_t nE = 20;
    // Sets the electric field range [V/cm] to be covered by the gas table. 
    const double Emin = 100., Emax = 100.E3;
    // Flag to request logarithmic spacing.
    constexpr bool useLog = true;

    // Gas pressure [Torr].
    const double gas_pressure = 1. * AtmosphericPressure;
    // Gas temperature [K]
    const double gas_temperature = 293.15;

    // .gas file name.
    const std::string gas_name = "Ar_75_iC4H10_25";
    // Saves the gas files here. 
    const std::string gas_folder = "./GasFiles/";  

    // Sets the gas composition, temperature and pressure.
    MediumMagboltz gas("Ar", 75., "iC4H10", 25.);
    gas.SetTemperature(gas_temperature);
    gas.SetPressure(gas_pressure);

    // Sets the points of evaluation.
    gas.SetFieldGrid(Emin, Emax, nE, useLog); 

    // Runs Magboltz to generate the gas table.
    gas.GenerateGasTable(nCollisions, true);

    // Merges the previous gas table with the current one if it exists.
    if(std::filesystem::exists(gas_folder + "/" + gas_name + ".gas")) {
        gas.MergeGasFile(gas_folder + "/" + gas_name + ".gas", true);
    }
    else {
        std::filesystem::create_directories(gas_folder);
    }

    // Saves the .gas file.
    gas.WriteGasFile(gas_folder + "/" + gas_name + ".gas");
}