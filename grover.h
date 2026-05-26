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

// Применение квантового преобразования Фурье (QFT) ко всем кубитам
// inverse = 0 -> прямое QFT, inverse = 1 -> обратное QFT
void apply_qft(SparseState *s, int inverse);

// Функция диффузии Гровера (отражение относительно среднего)
void apply_diffusion(SparseState *s);

int measure(SparseState *s, double *probability);

#endif