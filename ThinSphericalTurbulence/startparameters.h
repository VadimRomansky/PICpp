#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = true;
const bool turbulence = true;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::SPHERICAL;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW};
const Input input = Input::SMILEI;

enum Scale {LINEAR, LOG};
const Scale scale = Scale::LINEAR;

#endif