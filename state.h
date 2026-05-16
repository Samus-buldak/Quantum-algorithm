#ifndef STATE_H
#define STATE_H

#include <complex.h>

// Одна запись: номер базисного состояния и его комплексная амплитуда
typedef struct {
    int index;
    double complex amp;
} Entry;

// Разреженное состояние
typedef struct {
    int n_qubits;          // число кубитов
    int size;              // полное число состояний = 2^n_qubits
    Entry *entries;        // массив ненулевых записей
    int count;             // сколько сейчас записей
    double epsilon;        // порог отсечения (1e-10)
} SparseState;

// Создать состояние |0⟩ (одна запись с индексом 0 и амплитудой 1)
SparseState* create_state(int n_qubits);

// Освободить память
void free_state(SparseState *s);

// Получить амплитуду для заданного индекса (0 если нет)
double complex get_amp(SparseState *s, int index);

// Установить амплитуду для индекса (добавляет или удаляет запись, если амплитуда слишком мала)
void set_amp(SparseState *s, int index, double complex value);

// Вычислить сумму квадратов модулей (норму в квадрате)
double norm_sq(SparseState *s);

// Напечатать все ненулевые амплитуды (для отладки)
void print_state(SparseState *s);

#endif