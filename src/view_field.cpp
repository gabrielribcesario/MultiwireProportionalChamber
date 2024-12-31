#include <cstdlib>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ViewGeometry.hh>

using namespace Garfield;

int main(int argc, char* argv[]) {
    TApplication app("meshApp", &argc, argv);

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
    MediumMagboltz gas;

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
    
    ViewField fieldViewXZ(&mwpc);
    fieldViewXZ.SetArea(-anodeSpacing * 2., -acGap * 1.2, 
                         anodeSpacing * 2.,  acGap * 1.2);
    fieldViewXZ.SetPlaneXZ();
    fieldViewXZ.PlotContour();
    
    ViewField fieldViewYZ(&mwpc);
    fieldViewYZ.SetArea(-cathodeSpacing * 2., -acGap * 1.2, 
                         cathodeSpacing * 2.,  acGap * 1.2);
    fieldViewYZ.SetPlaneYZ();
    fieldViewYZ.PlotContour();

    ViewGeometry geoView(&mwpcGeo);
    geoView.SetArea(-cathodeSpacing * 2., -anodeSpacing * 2., -acGap * 1.2, 
                     cathodeSpacing * 2.,  anodeSpacing * 2.,  acGap * 1.2);
    geoView.Plot3d();

    app.Run(true);

    return 0;
}