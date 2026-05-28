#include "state.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

// Создание состояния |0⟩: для n<=20 сразу плотное, для n>20 разреженное
SparseState* create_state(int n_qubits) {
    SparseState *s = malloc(sizeof(SparseState));
    if (!s) return NULL;
    s->n_qubits = n_qubits;
    s->size = 1 << n_qubits;
    s->epsilon = 1e-10;
    s->is_mmap = 0;
    s->mmap_filename = NULL;

    if (n_qubits <= 20) {
        // Плотное состояние в оперативной памяти
        s->is_dense = 1;
        s->data.dense = calloc(s->size, sizeof(double complex));
        if (!s->data.dense) {
            free(s);
            return NULL;
        }
        s->data.dense[0] = 1.0;
        s->count = s->size;
    } else {
        // Для n > 20 используем mmap (дисковое хранилище)
        // Создаём временный файл с уникальным именем (можно фиксированное, но лучше mkstemp)
        char template[] = "state_mmap_XXXXXX";
        int fd = mkstemp(template);
        if (fd == -1) {
            perror("mkstemp");
            free(s);
            return NULL;
        }
        size_t file_size = (size_t)s->size * sizeof(double complex);
        if (ftruncate(fd, file_size) == -1) {
            perror("ftruncate");
            close(fd);
            free(s);
            return NULL;
        }
        void *ptr = mmap(NULL, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (ptr == MAP_FAILED) {
            perror("mmap");
            close(fd);
            free(s);
            return NULL;
        }
        close(fd); // файл можно закрыть, отображение остаётся
        // Заполняем нулями и устанавливаем |0⟩
        memset(ptr, 0, file_size);
        ((double complex*)ptr)[0] = 1.0;

        s->is_dense = 1;
        s->data.dense = ptr;
        s->count = s->size;
        s->is_mmap = 1;
        s->mmap_filename = strdup(template);
        if (!s->mmap_filename) {
            munmap(ptr, file_size);
            free(s);
            return NULL;
        }
    }
    return s;
}

void free_state(SparseState *s) {
    if (s) {
        if (s->is_dense) {
            if (s->is_mmap) {
                size_t size = (size_t)s->size * sizeof(double complex);
                munmap(s->data.dense, size);
                if (s->mmap_filename) {
                    remove(s->mmap_filename); // удаляем временный файл
                    free(s->mmap_filename);
                }
            } else {
                free(s->data.dense);
            }
        } else {
            free(s->data.entries);
        }
        free(s);
    }
}

// Линейный поиск
double complex get_amp(SparseState *s, int index) {
    if (s->is_dense) {                              // Для плотного представления                    
        if (index >= 0 && index < s->size)
            return s->data.dense[index];
        return 0.0;
    } else {                                        // Для разряженного
        for (int i = 0; i < s->count; i++)
            if (s->data.entries[i].index == index)
                return s->data.entries[i].amp;
        return 0.0;
    }
}

void set_amp(SparseState *s, int index, double complex value) {
    // Если амплитуда почти ноль – удаляем запись
    if (s->is_dense) {
        if (index >= 0 && index < s->size)
            s->data.dense[index] = value;
        return;
    }
    
    // Разреженный режим
    if (cabs(value) < s->epsilon) {
        // Удалить запись
        for (int i = 0; i < s->count; i++) {
            if (s->data.entries[i].index == index) {
                s->data.entries[i] = s->data.entries[s->count - 1];
                s->count--;
                return;
            }
        }
        return;
    }
    // Найти или добавить
    for (int i = 0; i < s->count; i++) {
        if (s->data.entries[i].index == index) {
            s->data.entries[i].amp = value;
            return;
        }
    }
    // Добавить новую
    if (s->count >= s->size) {
        // Эта ситуация не должна возникать для n>20, но на всякий случай переключимся в плотное?
        // Однако для n>20 плотное невозможно из-за памяти. Поэтому ошибка.
        fprintf(stderr, "Ошибка: разреженное состояние переполнено. n слишком велико для выполнения операций.\n");
        exit(1);
    }
    // Расширяем массив entries при необходимости
    if (s->count % 16 == 0 && s->count > 0) {
        Entry *new_entries = realloc(s->data.entries, (s->count + 16) * sizeof(Entry));
        if (!new_entries) { fprintf(stderr, "Ошибка памяти\n"); exit(1); }
        s->data.entries = new_entries;
    }
    s->data.entries[s->count].index = index;
    s->data.entries[s->count].amp = value;
    s->count++;
}

double norm_sq(SparseState *s) {
    double sum = 0.0;
    if (s->is_dense) {
        for (int i = 0; i < s->size; i++) {
            double m = cabs(s->data.dense[i]);
            sum += m * m;
        }
    } else {
        for (int i = 0; i < s->count; i++) {
            double m = cabs(s->data.entries[i].amp);
            sum += m * m;
        }
    }
    return sum;
}

void print_state(SparseState *s) {
    if (s->is_dense) {
        printf("Плотное состояние (все %d амплитуд):\n", s->size);
        for (int i = 0; i < s->size && i < 20; i++) { // печатаем только первые 20
            printf("|%d> : %.6f %+.6fi\n", i, creal(s->data.dense[i]), cimag(s->data.dense[i]));
        }
        if (s->size > 20) printf("... и ещё %d состояний\n", s->size - 20);
    } else {
        printf("Разреженное состояние: %d ненулевых записей\n", s->count);
        for (int i = 0; i < s->count; i++) {
            printf("|%d> : %.6f %+.6fi\n", s->data.entries[i].index,
                   creal(s->data.entries[i].amp), cimag(s->data.entries[i].amp));
        }
    }
}

double complex* sparse_to_dense_array(SparseState *s) {
    if (s->is_dense) {
        // Возвращаем копию плотного массива
        double complex *copy = malloc(s->size * sizeof(double complex));
        if (copy) memcpy(copy, s->data.dense, s->size * sizeof(double complex));
        return copy;
    } else {
        double complex *dense = calloc(s->size, sizeof(double complex));
        if (!dense) return NULL;
        for (int i = 0; i < s->count; i++)
            dense[s->data.entries[i].index] = s->data.entries[i].amp;
        return dense;
    }
}

void convert_to_dense(SparseState *s, double complex *dense) {
    // Принимает владение над массивом dense (не копирует)
    if (s->is_dense) {
        free(s->data.dense);
    } else {
        free(s->data.entries);
    }
    s->is_dense = 1;
    s->data.dense = dense;
    s->count = s->size;
}

int count_nonzero(SparseState *s) {
    if (s->is_dense) {
        int cnt = 0;
        for (int i = 0; i < s->size; i++) {
            if (cabs(s->data.dense[i]) >= s->epsilon) cnt++;
        }
        return cnt;
    } else {
        return s->count; // уже хранится
    }
}