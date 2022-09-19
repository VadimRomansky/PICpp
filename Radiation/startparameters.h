#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = false;
const bool turbulence = false;

const bool initialGridSearch = false;
const bool optimization = true; //false for debug only

enum Doppler {NO, INTEGER, DIFFERENTIAL};
const Doppler doppler = Doppler::DIFFERENTIAL;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::SPHERICAL;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW, COMBINED};
const Input input = Input::COMBINED;

enum Scale {LINEAR, LOG};
const Scale scale = Scale::LINEAR;

#endif