#ifndef CONSTANTS_H
#define CONSTANTS_H

const int randomParameter = 32*1024;
const double atomicUnionMass = 1.66053904E-24;
const double massProtonReal = 1.67262177E-24;
const double massAlphaReal = 6.644656E-24;
const double massElectron = 0.910938291E-27;
//const double massElectronFactor = massProtonReal/(10.0*massElectronReal);
const double massDeuteriumReal = 3.34449696893E-24;
const double massHelium3Real = 5.00823792874E-24;
const double massOxygenReal = 26.5601801672E-24;
const double massSiliconReal = 46.4567787264E-24;

const double massRelationSqrt = sqrt(100.0);
const double realMassRelationSqrt = sqrt(massProtonReal/massElectron);

const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double speed_of_light2 = speed_of_light * speed_of_light;
const double speed_of_light4 = speed_of_light2 * speed_of_light2;
const double electron_charge = 4.803529695E-10;
const double hplank = 6.626E-27;
const double pi = 4*atan2(1.0,1.0);
const double four_pi = 4*pi;

const double distance = 150*3*1.0E24;

const std::string outputfileName = "radiation.dat";
//const std::string fileNameP = "../../tristan-mp-pitp/Pe";
//const std::string fileNameF = "../../tristan-mp-pitp/Fe";
const std::string fileNameP = "Ee";
const std::string fileNameF = "Fs";
const std::string logFileName = "log.dat";


#endif