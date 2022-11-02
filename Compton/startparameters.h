#ifndef STARTPARAMETERS_H
#define STARTPARAMETERS_H

const bool parker = true;
const bool turbulence =false;

enum Geometry {FLAT_SIMPLE, FLAT, SPHERICAL};
const Geometry geometry = Geometry::SPHERICAL;

enum Input {TRISTAN, SMILEI, MAXWELL, POWERLAW, COMBINED, LONG_COMBINED};
const Input input = Input::POWERLAW;

enum Solver {DUBUS, UVAROV, THOMPSON};
const Solver solver = Solver::THOMPSON;

#endif