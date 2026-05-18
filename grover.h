#ifndef GROVER_H
#define GROVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "state.h"

// Функция применения гейта Адамара ко всем кубитам (H⊗n)
void apply_hadamard_all(SparseState *s);

#endif