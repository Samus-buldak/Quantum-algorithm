#include "grover.h"

// Функция сравнения для qsort (по индексу)
static int cmp_entry(const void *a, const void *b) {
    int ia = ((const Entry*)a)->index;
    int ib = ((const Entry*)b)->index;
    return (ia > ib) - (ia < ib);
}

void apply_hadamard_all(SparseState *s) {
    int n = s->n_qubits;

    // Если состояние плотное – используем быстрый O(N·n) алгоритм
    if (s->is_dense) {
        double factor = M_SQRT1_2;          // 1/√2
        int N = s->size;
        for (int qubit = 0; qubit < n; qubit++) {
            int step = 1 << qubit;
            for (int i = 0; i < N; i += 2 * step) {
                for (int j = 0; j < step; j++) {
                    int i0 = i + j;
                    int i1 = i + j + step;
                    double complex a0 = s->data.dense[i0];
                    double complex a1 = s->data.dense[i1];
                    s->data.dense[i0] = (a0 + a1) * factor;
                    s->data.dense[i1] = (a0 - a1) * factor;
                }
            }
        }
        return;
    }

    // ---------- Разреженный режим (оптимизированный) ----------
    double factor = M_SQRT1_2;

    for (int qubit = 0; qubit < n; qubit++) {
        int old_count = s->count;
        Entry *old = s->data.entries;

        // Временный массив: каждый старый индекс порождает 2 вклада
        int max_temp = old_count * 2;
        Entry *temp = malloc(max_temp * sizeof(Entry));
        if (!temp) {
            fprintf(stderr, "Ошибка памяти в apply_hadamard_all\n");
            return;
        }

        int tcnt = 0;
        for (int i = 0; i < old_count; i++) {
            int idx = old[i].index;
            double complex amp = old[i].amp;
            int bit = (idx >> qubit) & 1;
            int idx0 = idx;                     // бит = 0
            int idx1 = idx ^ (1 << qubit);      // бит = 1
            double complex c = amp * factor;

            if (bit == 0) {
                temp[tcnt].index = idx0; temp[tcnt].amp = c; tcnt++;
                temp[tcnt].index = idx1; temp[tcnt].amp = c; tcnt++;
            } else {
                temp[tcnt].index = idx0; temp[tcnt].amp = c; tcnt++;
                temp[tcnt].index = idx1; temp[tcnt].amp = -c; tcnt++;
            }
        }

        // Сортируем временный массив по индексу
        qsort(temp, tcnt, sizeof(Entry), cmp_entry);

        // Слияние записей с одинаковым индексом + отсечение малых амплитуд
        Entry *new_entries = malloc(tcnt * sizeof(Entry));
        int new_count = 0;
        for (int i = 0; i < tcnt; ) {
            int j = i;
            double complex sum = 0;
            while (j < tcnt && temp[j].index == temp[i].index) {
                sum += temp[j].amp;
                j++;
            }
            if (cabs(sum) >= s->epsilon) {
                new_entries[new_count].index = temp[i].index;
                new_entries[new_count].amp = sum;
                new_count++;
            }
            i = j;
        }

        free(temp);
        free(old);
        s->data.entries = new_entries;
        s->count = new_count;

        // (опционально) если число ненулевых записей превысило половину размера,
        // можно переключиться в плотный режим – для n ≤ 20
        if (!s->is_dense && s->count > s->size / 2 && s->size <= (1 << 20)) {
            double complex *dense = calloc(s->size, sizeof(double complex));
            if (dense) {
                for (int i = 0; i < s->count; i++)
                    dense[s->data.entries[i].index] = s->data.entries[i].amp;
                free(s->data.entries);
                s->is_dense = 1;
                s->data.dense = dense;
                s->count = s->size;
            }
        }
    }
}