#include <cstdlib>
#include <string>
#include <cmath>
#include <omp.h>

#include <TROOT.h>
#include <TObject.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TSystem.h>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/DriftLineRKF.hh>
#include <Garfield/TrackHeed.hh>

#include "CustomContainers.hh"
#include "Utility.hh"

// Relative probabilities per 100 disintegrations given the Fe-55 X-Ray emissions table.

#define KALPHA1_P 16.48  
#define KALPHA2_P 8.40
#define KBETA_P 3.38

// X-Ray energy [eV].

#define KALPHA1_E 5.89881E3
#define KALPHA2_E 5.88772E3
#define KBETA_E 6.49051E3

using namespace Garfield;

int main(int argc, char* argv[]) {
    // Loads the custom container library onto the ROOT system.
    gSystem->Load("libCustomContainers.so");

    const std::string resultsDir = InitializeDir("SimulationResults");

    // Creates a new .root results file. Compression algorithms: https://root.cern.ch/doc/master/Compression_8h_source.html
    // ZSTD is about as fast as LZ4 due to the ROOT IO API bottleneck, but has a much higher compression ratio.
    TFile resultsFile((resultsDir + "/SimulationOutput.root").c_str(), "create", "", 
                       ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);

    // Sets the gas and its composition.
    MediumMagboltz gas("Ar", 80.0, "CO2", 20.0);
    // Loads electron and ion transport parameters.
    gas.LoadGasFile("./GasFiles/Ar_80_CO2_20.gas");
    const std::string garfield_dir = std::getenv("GARFIELD_INSTALL");
    gas.LoadIonMobility(garfield_dir + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    // Enables Penning transfer. Uses a pre-implemented parameterization if no arguments are passed.
    gas.EnablePenningTransfer();

    // Anode-Cathode gap [cm];
    const double acGap = 0.5;
    // Wire length [cm].
    const double wireLength = 10.0;

    // Cathode wire spacing [cm].
    const double cathodeSpacing = 0.2;
    // Cathode wire diameter [cm].
    const double cathodeDiameter = 0.005;
    // Cathode wire potential [V].
    const double cathodeV = 0.0;

    // Anode wire spacing [cm].
    const double anodeSpacing = 0.4;
    // Anode wire diameter [cm].
    const double anodeDiameter = 0.002;
    // Anode wire potential [V].
    const double anodeV = 2100.0;

    MediumConductor metal;

    // Upper cathode wires plane, parallel to the x-axis.
    SolidWire cathodeWireU(0.0, 0.0, acGap, 
                           cathodeDiameter / 2.0, wireLength / 2.0, 
                           1.0, 0.0, 0.0);
    cathodeWireU.SetBoundaryPotential(cathodeV);
    //cathodeWireU.SetLabel("CathodeWireUpper");

    // Lower cathode wires plane, parallel to the x-axis.
    SolidWire cathodeWireL(0.0, 0.0, -acGap, 
                           cathodeDiameter / 2.0, wireLength / 2.0, 
                           1.0, 0.0, 0.0);
    cathodeWireL.SetBoundaryPotential(cathodeV);
    //cathodeWireL.SetLabel("CathodeWireLower");

    // Anode wires plane, parallel to the y-axis.
    SolidWire anodeWire(0.0, 0.0, 0.0, 
                        anodeDiameter / 2.0, wireLength / 2.0, 
                        0.0, 1.0, 0.0);
    anodeWire.SetBoundaryPotential(anodeV);
    //anodeWire.SetLabel("AnodeWire");

    GeometrySimple mwpcGeo;
    mwpcGeo.SetMedium(&gas);
    mwpcGeo.AddSolid(&cathodeWireU, &metal);
    mwpcGeo.AddSolid(&cathodeWireL, &metal);
    mwpcGeo.AddSolid(&anodeWire, &metal);

    // neBEM element size.
    const double elementSize = anodeDiameter < cathodeDiameter ? anodeDiameter / 10.0 : cathodeDiameter / 10.0;

    // Number of threads to use for parallel regions.
    const int n_jobs = omp_get_max_threads();

    ComponentNeBem3d mwpc;
    mwpc.SetGeometry(&mwpcGeo);
    mwpc.SetNumberOfThreads(n_jobs);
    mwpc.SetPeriodicityX(anodeSpacing);
    mwpc.SetPeriodicityY(cathodeSpacing);
    mwpc.SetTargetElementSize(elementSize);
    mwpc.UseLUInversion();
    mwpc.Initialise();

    // Creates the interface between the transport classes and the component.
    Sensor sensor;
    sensor.AddComponent(&mwpc);
    sensor.SetArea(-wireLength / 2.0, -wireLength / 2.0, -1.1 * (acGap + cathodeDiameter / 2.0),
                    wireLength / 2.0,  wireLength / 2.0,  1.1 * (acGap + cathodeDiameter / 2.0));

    // DriftLineRKF maximum step size.
    const double RKFStepSize = 0.01;
    // DriftLineRKF integration accuracy.
    const double RKFepsilon = 1.0E-8;

    std::printf("DriftLineRKF max step size [cm]: %#g\n", RKFStepSize);
    std::printf("DriftLineRKF integration accuracy: %#g\n", RKFepsilon);

    // Primary electron initial point, (x0, y0, z0, t0).
    CustomContainer::Position4D primaryInit;
    // Primary electron end point (x1, y1, z1, t1).
    CustomContainer::Position4D primaryEndP;
    // Primary electron status flag.
    int status = 0;
    // Arrival time spread [ns].
    double ats = 0.0;
    // Drift line gain.
    double gain = 0.0;
    // Drift line loss.
    double loss = 0.0;
    // # of electrons at the end of the drift line.
    double nEle = 0.0;
    // # of ions at the end of the ion tail.
    double nIon = 0.0;

    TTree driftData("DriftData", "Drift lines data");
    driftData.Branch("InitialPoint", &primaryInit);
    driftData.Branch("Endpoint", &primaryEndP);
    driftData.Branch("Status", &status);
    driftData.Branch("ArrivalTimeSpread", &ats);
    driftData.Branch("Gain", &gain);
    driftData.Branch("Loss", &loss);
    driftData.Branch("SizeElectrons", &nEle);
    driftData.Branch("SizeIons", &nIon);

    // Class responsible for the primary ionizations calculations.
    TrackHeed track(&sensor);

    const size_t nEvents = 1000;

    TH1D trackHist("PhotonConversion", "X-ray Energy Spectrum", 100, 0, 250);

    for (unsigned int i = 0; i < nEvents; ++i) {  
        // Photon initial point, (x0, y0, z0, t0).
        CustomContainer::Position4D photonInit;
        // Photon initial direction, (dx, dy, dz).
        CustomContainer::Direction3D photonDir;
        // Photon energy [eV].
        double photonEnergy = 0.0;
        // Number of electron-ion pairs produced by the photon.
        int nPel = 0;  

        // Discards photons that do not interact with the detector.
        while (!nPel) {
            photonInit.SetValue(Uniform(-0.05, 0.05), Uniform(-0.05, 0.05), 1.1 * (acGap + cathodeDiameter / 2.0), 0.0);
            photonDir.SetValue(Uniform(-0.05, 0.05), Uniform(-0.05, 0.05), -1.0);

            double prob = Uniform(0.0, KALPHA1_P + KALPHA2_P + KBETA_P);
            photonEnergy = prob < KBETA_P ? KBETA_E : prob < KBETA_P + KALPHA1_P ? KALPHA1_E : KALPHA2_E; 

            track.TransportPhoton(photonInit.GetX(), photonInit.GetY(), photonInit.GetZ(), photonInit.GetT(),
                                  photonEnergy, 
                                  photonDir.GetDX(), photonDir.GetDY(), photonDir.GetDZ(),
                                  nPel);
        }

        trackHist.Fill(nPel);

        std::printf("Event #%d\n", i + 1);
        std::printf("|  Photon energy [eV]: %#g\n", photonEnergy);
        std::printf("|  Photon conversion: %i electron-ion pairs\n", nPel);
        std::printf("|  (x0, y0, z0): (%.6f, %.6f, %.6f)\n", photonInit.GetX(), photonInit.GetY(), photonInit.GetZ());
        std::printf("|  (dx, dy, dz): (%.6f, %.6f, %.6f)\n", photonDir.GetDX(), photonDir.GetDY(), photonDir.GetDZ());

        // Loops over the electrons produced by the photon's trajectory.
        for (const auto& cluster : track.GetClusters()) {
            #pragma omp parallel for schedule(dynamic, 1) num_threads(n_jobs)
            for (const auto& electron : cluster.electrons) {
                // Creates the class responsible for tracking secondary ionizations and connects it to a sensor.
                DriftLineRKF driftRKF(&sensor);
                driftRKF.SetMaximumStepSize(RKFStepSize);
                driftRKF.SetIntegrationAccuracy(RKFepsilon);

                double x0 = electron.x, y0 = electron.y, z0 = electron.z, t0 = electron.t;
                driftRKF.DriftElectron(x0, y0, z0, t0);

                int status_i = 0;
                double x1 = 0.0, y1 = 0.0, z1 = 0.0, t1 = 0.0;
                driftRKF.GetEndPoint(x1, y1, z1, t1, status_i);

                double nEle_i = 0.0, nIon_i = 0.0;
                driftRKF.GetAvalancheSize(nEle_i, nIon_i);

                // Integration outside the critical region for better leveraging the parallel loop.
                double ats_i = driftRKF.GetArrivalTimeSpread(), gain_i = driftRKF.GetGain(), loss_i = driftRKF.GetLoss();

                #pragma omp critical 
                {
                    primaryInit.SetValue(x0, y0, z0, t0);
                    primaryEndP.SetValue(x1, y1, z1, t1);
                    status = status_i;
                    ats = ats_i;
                    gain = gain_i;
                    loss = loss_i;
                    nEle = nEle_i;
                    nIon = nIon_i;

                    driftData.Fill();
                }
            }
        }
        driftData.Write(nullptr, TObject::kWriteDelete);
        trackHist.Write(nullptr, TObject::kWriteDelete);
    }
    resultsFile.Close();

    std::printf("Done.\n");

    return EXIT_SUCCESS;
}