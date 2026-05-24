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

// Вспомогательная функция: БПФ (бабочка) для плотного массива
// sign = 1 для ядра exp(+2πi·...), sign = -1 для exp(-2πi·...)
static void fft_core(double complex *state, int n, int sign) {
    int N = 1 << n;
    // Битово-инверсная перестановка
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int b = 0; b < n; b++) {
            if (i & (1 << b))
                rev |= (1 << (n - 1 - b));
        }
        if (i < rev) {
            double complex tmp = state[i];
            state[i] = state[rev];
            state[rev] = tmp;
        }
    }

    // Итеративные бабочки
    for (int len = 2; len <= N; len <<= 1) {
        double angle = sign * 2.0 * M_PI / len;
        double complex wlen = cos(angle) + I * sin(angle);
        for (int i = 0; i < N; i += len) {
            double complex w = 1.0;
            for (int j = 0; j < len / 2; j++) {
                double complex u = state[i + j];
                double complex v = state[i + j + len / 2] * w;
                state[i + j] = u + v;
                state[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// Прямое или обратное квантовое преобразование Фурье
void apply_qft(SparseState *s, int inverse) {
    if (!s) return;

    // Для n > 20 QFT не реализовано в разреженном виде (станет плотным)
    if (!s->is_dense) {
        // Если разреженное состояние, но число кубитов ≤20 – преобразуем в плотное
        if (s->n_qubits <= 20 && s->size <= (1 << 20)) {
            double complex *dense = sparse_to_dense_array(s);
            if (!dense) {
                fprintf(stderr, "Ошибка памяти при преобразовании в плотное состояние\n");
                return;
            }
            convert_to_dense(s, dense);
        } else {
            fprintf(stderr, "QFT поддерживается только для плотного состояния (n ≤ 20)\n");
            return;
        }
    }

    // Теперь состояние плотное
    double complex *state = s->data.dense;
    int n = s->n_qubits;
    int N = s->size;

    // Выбор знака: прямое QFT использует exp(+2πi·xy/N), обратное – exp(-2πi·xy/N)
    int sign = inverse ? -1 : 1;
    fft_core(state, n, sign);

    // Нормировка: обе формулы содержат множитель 1/√N
    double inv_sqrtN = 1.0 / sqrt(N);
    for (int i = 0; i < N; i++) {
        state[i] *= inv_sqrtN;
    }
}

void apply_diffusion(SparseState *s) {
    if (!s) return;

    // Для n > 20 разреженное состояние не поддерживает эффективную диффузию
    if (!s->is_dense) {
        if (s->n_qubits <= 20) {
            // Преобразуем в плотное (станет размер 2^n, для n≤20 это OK)
            double complex *dense = sparse_to_dense_array(s);
            if (!dense) {
                fprintf(stderr, "Ошибка памяти в apply_diffusion\n");
                return;
            }
            convert_to_dense(s, dense);
        } else {
            fprintf(stderr, "Диффузия не поддерживается для разреженного состояния с n > 20\n");
            return;
        }
    }

    // Теперь состояние плотное
    int N = s->size;
    double complex *a = s->data.dense;

    // Вычисляем среднюю амплитуду
    double complex sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += a[i];
    }
    double complex avg = sum / N;

    // Отражение относительно среднего
    for (int i = 0; i < N; i++) {
        a[i] = 2.0 * avg - a[i];
    }
}