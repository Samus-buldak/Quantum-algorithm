#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "state.h"

// Функция применения гейта Адамара ко всем кубитам (H⊗n)
void apply_hadamard_all(SparseState *s) {
    int N = s->size;
    // Временные полные векторы (пока что используем полное представление для простоты)
    double complex *old = (double complex*) malloc(N * sizeof(double complex));
    double complex *new = (double complex*) malloc(N * sizeof(double complex));
    if (!old || !new) {
        printf("Ошибка выделения памяти\n");
        free(old); free(new);
        return;
    }
    // Инициализируем old нулями
    for (int i = 0; i < N; i++) old[i] = 0.0;
    // Копируем ненулевые амплитуды
    for (int i = 0; i < s->count; i++) {
        int idx = s->entries[i].index;
        old[idx] = s->entries[i].amp;
    }
    double factor = 1.0 / sqrt((double)N);
    // Вычисляем new[x] = (1/√N) * Σ_y (-1)^{x·y} old[y]
    for (int x = 0; x < N; x++) {
        double complex sum = 0.0;
        for (int y = 0; y < N; y++) {
            // Скалярное произведение битов x и y по модулю 2
            int bits = x & y;
            int dot = 0;
            while (bits) {
                dot ^= (bits & 1);
                bits >>= 1;
            }
            double sign = (dot == 0) ? 1.0 : -1.0;
            sum += sign * old[y];
        }
        new[x] = sum * factor;
    }
    // Перезаписываем разреженное состояние
    s->count = 0;
    for (int x = 0; x < N; x++) {
        if (cabs(new[x]) >= s->epsilon) {
            set_amp(s, x, new[x]);
        }
    }
    free(old);
    free(new);
}

int main() {
    int n = 4;  // 3 кубита → 8 состояний
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