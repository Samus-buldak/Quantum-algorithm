#include "state.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Создать состояние |0⟩
SparseState* create_state(int n_qubits) {
    SparseState *s = (SparseState*) malloc(sizeof(SparseState));
    if (!s) return NULL;
    s->n_qubits = n_qubits;
    s->size = 1 << n_qubits;
    s->epsilon = 1e-10;
    // Начальная вместимость – 16 записей, но для простоты выделим сразу максимально возможное количество
    // (на небольших n (≤20) это допустимо, память ~ 2^20 * 16 байт ≈ 16 МБ)
    // Если хотите экономию, можно выделять постепенно, но для простоты оставим так.
    s->entries = (Entry*) malloc(s->size * sizeof(Entry));
    if (!s->entries) {
        free(s);
        return NULL;
    }
    s->count = 0;
    // Устанавливаем |0⟩
    set_amp(s, 0, 1.0);
    return s;
}

void free_state(SparseState *s) {
    if (s) {
        free(s->entries);
        free(s);
    }
}

// Линейный поиск (для простоты, для n≤20 скорость достаточная)
double complex get_amp(SparseState *s, int index) {
    for (int i = 0; i < s->count; i++) {
        if (s->entries[i].index == index)
            return s->entries[i].amp;
    }
    return 0.0;
}

void set_amp(SparseState *s, int index, double complex value) {
    // Если амплитуда почти ноль – удаляем запись
    if (cabs(value) < s->epsilon) {
        for (int i = 0; i < s->count; i++) {
            if (s->entries[i].index == index) {
                // замена последним элементом
                s->entries[i] = s->entries[s->count - 1];
                s->count--;
                return;
            }
        }
        return;
    }
    // Ищем существующую запись
    for (int i = 0; i < s->count; i++) {
        if (s->entries[i].index == index) {
            s->entries[i].amp = value;
            return;
        }
    }
    // Добавляем новую
    if (s->count >= s->size) {
        // Такое может случиться только если все амплитуды ненулевые – тогда разреженность теряется.
        // Для простоты просто выходим (в реальном симуляторе нужно увеличить массив)
        return;
    }
    s->entries[s->count].index = index;
    s->entries[s->count].amp = value;
    s->count++;
}

double norm_sq(SparseState *s) {
    double sum = 0.0;
    for (int i = 0; i < s->count; i++) {
        double m = cabs(s->entries[i].amp);
        sum += m * m;
    }
    return sum;
}

void print_state(SparseState *s) {
    printf("Ненулевых записей: %d\n", s->count);
    for (int i = 0; i < s->count; i++) {
        printf("|%d> : %.6f %+.6fi\n",
               s->entries[i].index,
               creal(s->entries[i].amp),
               cimag(s->entries[i].amp));
    }
}