#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = true;
const bool turbulence =false;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::SPHERICAL;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW};
const Input input = Input::SMILEI;

enum Solver {DUBUS, UVAROV};
const Solver solver = Solver::UVAROV;

#endif