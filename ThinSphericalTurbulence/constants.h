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
const double dzeta3 = 1.202056903;

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

const double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / (massElectron * speed_of_light2);
const double absorpCoef = 40.0*electron_charge /(9*pow(2, 4.0/3.0));
const double absorpCoef2 = 2.0 * electron_charge / 3.0 ;
const double criticalNuCoef = 3 * electron_charge / (4 * pi * massElectron * massElectron * massElectron * speed_of_light * speed_of_light4);

const std::string outputfileName = "radiation.dat";
//const std::string fileNameP = "../../tristan-mp-pitp/Pe";
//const std::string fileNameF = "../../tristan-mp-pitp/Fe";
const std::string fileNameP = "Ee";
const std::string fileNameF = "Fs";
const std::string logFileName = "log.dat";
const std::string BFileName = "B.dat";

const int Niterations = 20;
const int Nmonth = 6;
const int Nmontecarlo = 100000;

const double size[Nmonth] = {3.4E16, 6.8E16, 1.06E17, 2.4E17, 4.6E17, 14E17};
//const double size[Nmonth] = {3.4E16, 6.8E16, 1.06E17, 2.4E17};


const double times[Nmonth] = {0, 2760000, 5270400, 10700000, 19000000, 56100000};
//const double times[Nmonth] = {0, 2760000, 5270400, 10700000};

const double distance = 40*3.08*1.0E24;

const int Nrho = 10;
const int Nphi = 10;
const int Nz = 10;

const int Nk = 10;

//const double rmax = 3.4E16;

//const double dr = rmax/Nrho;
//const double dz = dr;

const double dphi = (2*pi)/Nphi;

const int Ntheta = 10;

const double dtheta = (pi/2)/Ntheta;
const double thetaValue[Ntheta] = {dtheta/2, 3*dtheta/2, 5*dtheta/2, 7*dtheta/2, 9*dtheta/2, 11*dtheta/2, 13*dtheta/2, 15*dtheta/2, 17*dtheta/2, 19*dtheta/2};
const double sinThetaValue[Ntheta] = {sin(dtheta/2), sin(3*dtheta/2), sin(5*dtheta/2), sin(7*dtheta/2), sin(9*dtheta/2), sin(11*dtheta/2), sin(13*dtheta/2), sin(15*dtheta/2), sin(17*dtheta/2), sin(19*dtheta/2)};
const double cosThetaValue[Ntheta] = {cos(dtheta/2), cos(3*dtheta/2), cos(5*dtheta/2), cos(7*dtheta/2), cos(9*dtheta/2), cos(11*dtheta/2), cos(13*dtheta/2), cos(15*dtheta/2), cos(17*dtheta/2), cos(19*dtheta/2)};


const int Napprox = 55;

///// x
const double UvarovX[Napprox] = {1.0E-9, 5.0E-9, 1.0E-8, 5.0E-8, 1.0E-7, 5.0E-7, 1.0E-6, 5.0E-6,
	0.00001, 0.00005, 0.0001 ,0.0005, 0.001, 0.005, 0.01, 0.025, 0.050, 0.075, 0.10, 0.15,
	0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.2, 1.4, 1.6,
	1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
	25.0, 30.0, 40.0, 50.0
};

///      K1/3(x)
const double McDonaldValue1_3[Napprox] = { 1687, 986.9, 783.3, 458.08, 363.58, 212.6, 168.7, 98.66,
	78.29, 45.74, 36.28, 21.134, 16.715, 9.5937, 7.486, 5.3015, 3.991, 3.3296, 2.8998, 2.3428, 1.9793, 1.71449, 1.509, 1.20576, 0.9890, 0.82509, 0.69653, 0.59318, 0.50859, 0.438,
	0.3298, 0.25132, 0.19324, 0.149649, 0.1165, 0.06354, 0.0353059, 0.01987, 0.011299, 0.006472, 0.003728, 0.0012547, 0.000427968, 0.000147, 0.00005118, 0.00001787, 0.00000221, 2.7719E-7, 3.511E-8, 4.4822E-9, 5.7568E-10,
	3.4717E-12, 2.13636E-14, 8.40438E-19, 3.4139E-23
};
///      K2/3(x)
const double McDonaldValue2_3[Napprox] = { 1.0747E6, 3.675E5, 2.315E5, 79189, 49886, 17060, 10747, 3675,
	2315, 791.89, 498.85, 170.6, 107.46, 36.72, 23.098, 12.468, 7.7619, 5.84349, 4.7529, 3.51287, 2.80179, 2.3289, 1.9866, 1.5171, 1.2059, 0.98283, 0.81478, 0.68387, 0.57938, 0.494,
	0.36618, 0.2756, 0.2099, 0.16132, 0.1248, 0.06725, 0.037057, 0.02073, 0.01173, 0.00669, 0.00384, 0.001287, 0.0004376, 0.0001503, 0.00005208, 0.00001816, 0.00000224, 2.80406E-7, 3.5469E-8, 4.5227E-9, 5.8038E-10,
	3.49449E-12, 2.14807E-14, 8.43904E-19, 3.4252E-23
};
///      K4/3(x)
const double McDonaldValue4_3[Napprox] = { 1.125E12, 1.3159E11, 5.222E10, 6.107E9, 2.4239E9, 2.835E8, 1.125E8, 1.315E7,
	5.222E6, 6.107E5, 2.4239E5, 28350, 11250.8, 1315.88, 522.179, 153.84, 60.97, 35.44, 24.085, 13.925, 9.39959, 6.90088, 5.34, 3.5267, 2.5246, 1.8996, 1.478, 1.178188, 0.9561, 0.78676,
	0.5494, 0.3953, 0.29044, 0.2167, 0.16368, 0.08419, 0.0449028, 0.0245, 0.01361, 0.00765, 0.00434, 0.0014269, 0.000478, 0.0001626, 0.00005587, 0.0000193, 0.000002363, 2.93606E-7, 3.693E-8, 4.6888E-9, 5.9957E-10,
	3.58707E-12, 2.19555E-14, 8.579119E-19, 3.47072E-23
};

///      K5/3(x)
const double McDonaldValue5_3[Napprox] = { 1.43E15, 9.8E13, 3.09E13, 2.11E12, 6.65E11, 4.55E10, 1.43E10, 9.8E8,
	3.08E8, 2.11E7, 6.65E6, 4.55E5, 1.43E4, 9802, 3087, 670, 211, 107, 66.3, 33.6, 20.7, 14.13, 10.338, 6.2628, 4.2048, 3.00916, 2.2484, 1.73296, 1.3669, 1.0977,
	0.73675, 0.513823, 0.36818, 0.269146, 0.1997, 0.0994, 0.051775, 0.02777, 0.0152, 0.008455, 0.004754, 0.0015408, 0.000511, 0.00017249, 0.00005889, 0.0000203, 0.00000246, 3.039E-7, 3.8068E-8, 4.817E-9, 6.1437E-10,
	3.658E-12, 2.2318E-14, 8.68568E-19, 3.50526E-23
};
////     x*int from x to inf K5/3(t)dt
const double UvarovValue[Napprox] = {0.0021495, 0.00367, 0.00463, 0.00792, 0.00997, 0.017, 0.0215, 0.0367,
	0.0461, 0.0791, 0.0995, 0.169, 0.213, 0.358, 0.445, 0.583, 0.702, 0.772, 0.818, 0.874,0.904, 0.917, 0.918, 0.901, 0.872, 0.832, 0.788,
	0.742, 0.694, 0.655, 0.566, 0.486, 0.414, 0.354, 0.301, 0.200, 0.130, 0.0845, 0.0541, 0.0339, 0.0214, 0.0085, 0.0033, 0.0013, 0.00050, 0.00019,
	0.00002822, 0.00000409, 5.89E-7, 8.42E-8, 1.19E-8,
	8.9564E-11, 6.58079E-13, 3.42988E-17, 1.73478E-21
};

const int Npsi = 6;
const int Nalpha = 6;
const double dalpha = (pi/2)/Nalpha;
const double alphaValue[Nalpha] = {dalpha/2, 3*dalpha/2, 5*dalpha/2, 7*dalpha/2, 9*dalpha/2, 11*dalpha/2};
const double sinAlphaValue[Nalpha] = {sin(dalpha/2), sin(3*dalpha/2), sin(5*dalpha/2), sin(7*dalpha/2), sin(9*dalpha/2), sin(11*dalpha/2)};
const double cosAlphaValue[Nalpha] = {cos(dalpha/2), cos(3*dalpha/2), cos(5*dalpha/2), cos(7*dalpha/2), cos(9*dalpha/2), cos(11*dalpha/2)};

const double decx[3]= {0.325, 0.61, 1.28};
const double decy[3] = {6.1, 1.9, 0.9};

const double octx[3]= {0.325, 0.61, 1.28};
const double octy[3] = {12.0, 6.3, 4.4};

const double augx[6] = {0.332, 0.617, 1.43, 4.86, 8.46, 84.6};
const double augy[6] = {3.3, 7.9, 8.68, 2.47, 1.084, 0.1084};

const double augmaxx = 0.8;
const double augmaxy = 11;

const double junx[4] = {0.617, 1.43, 4.86, 8.46};
const double juny[4] = {2.98, 12.3, 5.79, 3.15};

const double junmaxx = 1.65;
const double junmaxy = 13.2;

const double mayx[3]= {1.43, 4.86, 8.46};
const double mayy[3] = {4.93, 12.2, 6.82};

const double maymaxx = 2.96;
const double maymaxy = 15.2;

const double aprx[4] = {1.43, 4.86, 8.46, 22.5};
const double apry[4] ={1.3, 12.86, 17.57, 5.2};

const double aprmaxx = 6.50;
const double aprmaxy = 19.3;

const double minB = 0.0001;
const double maxB = 10;
const double minN = 0.01;
const double maxN = 400;
const double minFraction = 0.001;
const double maxFraction = 0.2;
const double maxSigma = 10.0;
const double minEta = 0;
const double maxEta = 0.8*pi/2;
const double maxV = 0.9*speed_of_light;
const double minV = 0.6*speed_of_light;
const double maxR = 30*24*3600*3E10;
const double minR = 1E15;

#endif