#include <cstdlib>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/DriftLineRKF.hh>
#include <Garfield/TrackHeed.hh>
#include <Garfield/ViewDrift.hh>
#include <Garfield/Random.hh>

#include "CustomContainers.hpp"

using namespace Garfield;

// Samples a value from a uniform distribution X ~ U(a,b).
double Uniform(double a, double b) { return (a - b) * RndmUniform() - a; }

int main(int argc, char* argv[]) {
    // Loads the custom container library onto the ROOT system.
    gSystem->Load("libCustomContainers.so");

    TApplication app("mwpcApp", &argc, argv);

    // Sets the gas and its composition.
    MediumMagboltz gas("Ar", 80., "CO2", 20.);
    // Loads electron and ion transport parameters.
    gas.LoadGasFile("./GasFiles/Ar_80_CO2_20.gas");
    const std::string garfield_dir = std::getenv("GARFIELD_INSTALL");
    gas.LoadIonMobility(garfield_dir + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    // Enables Penning transfer. Uses a pre-implemented parameterization if no arguments are passed.
    gas.EnablePenningTransfer();
    // Prepares the electron scattering rates table for the current gas mixture and density.
    gas.Initialise(true); 

    // Anode-Cathode gap [cm];
    const double acGap = 0.5;
    // Wire length [cm].
    const double wireLength = 10.;

    // Cathode wire spacing [cm].
    const double cathodeSpacing = 0.2;
    // Cathode wire diameter [cm].
    const double cathodeDiameter = 0.005;
    // Cathode wire potential [V].
    const double cathodeV = 0.;

    // Anode wire spacing [cm].
    const double anodeSpacing = 0.4;
    // Anode wire diameter [cm].
    const double anodeDiameter = 0.002;
    // Anode wire potential [V].
    const double anodeV = -2100.;

    MediumConductor metal;

    // Upper cathode wires plane, parallel to the z-axis.
    SolidWire cathodeWireU(0., 0., acGap, 
                           cathodeDiameter / 2., wireLength / 2., 
                           1., 0., 0.);
    cathodeWireU.SetBoundaryPotential(cathodeV);
    cathodeWireU.SetLabel("CathodeWireUpper");

    // Lower cathode wires plane, parallel to the z-axis.
    SolidWire cathodeWireL(0., 0., -acGap, 
                           cathodeDiameter / 2., wireLength / 2., 
                           1., 0., 0.);
    cathodeWireL.SetBoundaryPotential(cathodeV);
    cathodeWireL.SetLabel("CathodeWireLower");

    // Anode wires plane, parallel to the x-axis.
    SolidWire anodeWire(0., 0., 0., 
                        anodeDiameter / 2., wireLength / 2., 
                        0., 1., 0.);
    anodeWire.SetBoundaryPotential(anodeV);
    anodeWire.SetLabel("AnodeWire");

    GeometrySimple mwpcGeo;
    mwpcGeo.SetMedium(&gas);
    mwpcGeo.AddSolid(&cathodeWireU, &metal);
    mwpcGeo.AddSolid(&cathodeWireL, &metal);
    mwpcGeo.AddSolid(&anodeWire, &metal);
    
    // neBEM element size.
    const double elementSize = anodeDiameter < cathodeDiameter ? anodeDiameter : cathodeDiameter;

    ComponentNeBem3d mwpc;
    mwpc.SetGeometry(&mwpcGeo);
    mwpc.SetPeriodicityX(anodeSpacing);
    mwpc.SetPeriodicityY(cathodeSpacing);
    mwpc.SetTargetElementSize(elementSize);
    mwpc.UseLUInversion();
    mwpc.SetNumberOfThreads(8);
    mwpc.Initialise();

    // Creates the interface between the transport classes and the component.
    Sensor sensor;
    sensor.AddComponent(&mwpc);
    sensor.AddElectrode(&mwpc, "CathodeWireUpper");
    sensor.AddElectrode(&mwpc, "CathodeWireLower");
    sensor.AddElectrode(&mwpc, "AnodeWire");
    sensor.SetArea(-wireLength / 2., -wireLength / 2., -1.1 * (acGap + cathodeDiameter / 2.),
                    wireLength / 2.,  wireLength / 2.,  1.1 * (acGap + cathodeDiameter / 2.));

    // The ionizing particle traversing the gas chamber.
    const std::string trackParticle = "muon";
    // Particle momentum [eV/c].
    const double trackMomentum = 1.E9;
    // Particle initial point, (x0, y0, z0, t0).
    CustomContainer::Position4D trackInit(0., 0., 1.1 * (acGap + cathodeDiameter / 2.), 0.);
    // Particle initial direction, (dx, dy, dz).
    CustomContainer::Direction3D trackDir(Uniform(-1., 1.), Uniform(-1., 1.), -1.);

    // Class responsible for the primary ionizations calculations.
    TrackHeed track(&sensor);
    track.DisableDeltaElectronTransport();
    track.SetParticle(trackParticle);
    track.SetMomentum(trackMomentum);
    track.NewTrack(trackInit.GetX(), trackInit.GetY(), trackInit.GetZ(), trackInit.GetT(),
                   trackDir.GetDX(), trackDir.GetDY(), trackDir.GetDZ());

    std::printf("Particle: %s\n", trackParticle.c_str());
    std::printf("Momentum [eV/c]: %#g\n", trackMomentum);
    std::printf("(x0, y0, z0): (%.6f, %.6f, %.6f)\n", trackInit.GetX(), trackInit.GetY(), trackInit.GetZ());
    std::printf("(dx, dy, dz): (%.6f, %.6f, %.6f)\n", trackDir.GetDX(), trackDir.GetDY(), trackDir.GetDZ());
    std::printf("%ld electron clusters generated.\n", track.GetClusters().size());

    ViewDrift driftView;

    // DriftLineRKF step size.
    const double RKFStepSize = 1.E-3;
    // DriftLineRKF integration accuracy.
    const double RKFepsilon = 1.E-7;
    
    // Creates the class responsible for tracking secondary ionizations and connects it to a sensor.
    DriftLineRKF driftRKF(&sensor);
    driftRKF.EnableIonTail(true);
    driftRKF.EnableSignalCalculation(true);
    driftRKF.SetIntegrationAccuracy(RKFepsilon);
    driftRKF.SetMaximumStepSize(RKFStepSize);
    driftRKF.EnablePlotting(&driftView);

    std::printf("DriftLineRKF max step size [cm]: %.6f\n", RKFStepSize);
    std::printf("DriftLineRKF integration accuracy: %.6f\n", RKFepsilon);

    // Creates a new .root results file. Compression algorithms: https://root.cern.ch/doc/master/Compression_8h_source.html
    // ZSTD is about as fast as LZ4 due to the ROOT IO API bottleneck, but has a much higher compression ratio.
    TFile eventResults("./BetaParticleEvent.root", "recreate", "", 
                       ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);

    // Primary electrons data tree
    TTree primaryTree("PrimaryElectrons", "Primary Ionization Electrons Data");
    CustomContainer::Position4D primaryInit, primaryEndP;
    double gain, loss, ats;
    double nEle, nIon;
    int status;
    primaryTree.Branch("InitialPoint", &primaryInit, "x0/D:y0:z0:t0");
    primaryTree.Branch("EndPoint", &primaryEndP, "x1/D:y1:z1:t1");
    primaryTree.Branch("Status", &status);
    primaryTree.Branch("ArrivalTimeSpread", &ats);
    primaryTree.Branch("Gain", &gain);
    primaryTree.Branch("Loss", &loss);
    primaryTree.Branch("SizeElectrons", &nEle);
    primaryTree.Branch("SizeIons", &nIon);

    // Loops over the electron clusters along the particle's trajectory.
    for (const auto& cluster : track.GetClusters()) {
        // Loops over the electrons in the cluster.
        for (const auto& cstrElectron : cluster.electrons) {
            // Instantiates the tracking and avalanching of the electron.
            driftRKF.DriftElectron(cstrElectron.x, cstrElectron.y, cstrElectron.z, cstrElectron.t);

            double xx, yy, zz, tt;
            driftRKF.GetEndPoint(xx, yy, zz, tt, status);

            // Computes the gain, the loss and the arrival time spread outside of the critical section for faster execution.
            gain = driftRKF.GetGain();
            loss = driftRKF.GetLoss();
            ats = driftRKF.GetArrivalTimeSpread();
            driftRKF.GetAvalancheSize(nEle, nIon);
            primaryInit.setValue(cstrElectron.x, cstrElectron.y, cstrElectron.z, cstrElectron.t);
            primaryEndP.setValue(xx, yy, zz, tt);

            // Writes the collected data to their respective branches.
            primaryTree.Fill();
        }
    }
    primaryTree.Write(nullptr, TObject::kWriteDelete);
    eventResults.Close();

    TCanvas cSignalA("signal", "", 600, 600);
    sensor.PlotSignal("AnodeWire", &cSignalA);

    TCanvas cSignalCU("signal", "", 600, 600);
    sensor.PlotSignal("CathodeWireUpper", &cSignalCU);

    TCanvas cSignalCL("signal", "", 600, 600);
    sensor.PlotSignal("CathodeWireLower", &cSignalCL);

    std::printf("Done.\n");

    app.Run();

    return 0;
}