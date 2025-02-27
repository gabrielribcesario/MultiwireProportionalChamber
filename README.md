# Multiwire Proportional Chamber

This repository is meant to document the use of Garfield++ and PyROOT. For this project I chose one of NASA's experiments in trajectory measurement. 

Garfield's periodic boundary condition does not allow for the measurement of the induced signal in each individual wire, which means that the mesh I created using neBEM cannot be used to measure a particle's trajectory to the same degree of accuracy as the real chamber in a timely manner. It can however be used to measure the chamber's gain and, much like in the Garfield++ Heed example, the Fe-55 X-Ray energy spectrum.

/lib/ libraries might require recompilation depending on the ROOT version.