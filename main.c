#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "state.h"
#include "grover.h"

int main() {
    int n = 3;  // 3 кубита → 8 состояний
    printf("=== Симуляция квантового регистра с %d кубитами ===\n", n);
    SparseState *state = create_state(n);
    if (!state) {
        printf("Ошибка создания состояния\n");
        return 1;
    }
    printf("\nНачальное состояние |0>:\n");
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));
    
    printf("\nПрименяем гейт Адамара H⊗%d:\n", n);
    apply_hadamard_all(state);
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));
    
    free_state(state);
    return 0;
}