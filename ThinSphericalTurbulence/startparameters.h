#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = false;
const bool turbulence = false;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::FLAT;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW};
const Input input = Input::SMILEI;

enum Scale {LINEAR, LOG};
const Scale scale = Scale::LINEAR;

#endif