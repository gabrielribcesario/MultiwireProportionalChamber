#ifdef __CLING__

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ namespace CustomContainer;

#pragma link C++ class CustomContainer::Position4D+;
#pragma link C++ class CustomContainer::Direction3D+;

#pragma link C++ class vector< CustomContainer::Position4D >+;
#pragma link C++ class vector< CustomContainer::Direction3D >+;

#pragma link C++ class ROOT::VecOps::RVec< CustomContainer::Position4D >+;
#pragma link C++ class ROOT::VecOps::RVec< CustomContainer::Direction3D >+;

#endif