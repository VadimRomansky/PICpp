#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = false;
const bool turbulence = false;

enum Doppler { NO, INTEGER, DIFFERENTIAL };
const Doppler doppler = Doppler::INTEGER;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::SPHERICAL;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW, COMBINED, MONTECARLO, MONTECARLO_CUT, MONOENERGETIC};
const Input input = Input::COMBINED;

enum Scale {LINEAR, LOG};
const Scale scale = Scale::LINEAR;

#endif