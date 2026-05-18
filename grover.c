#include "grover.h"

// Функция применения гейта Адамара ко всем кубитам (H⊗n)
static void add_amp(SparseState *s, int index, double complex value) {
    double complex old = get_amp(s, index);
    set_amp(s, index, old + value);
}

void apply_hadamard_all(SparseState *s) {
    int n = s->n_qubits;
    double factor = 1.0 / sqrt(2.0);
    
    // Проходим по каждому кубиту
    for (int qubit = 0; qubit < n; qubit++) {
        // Копируем текущие ненулевые записи во временный массив
        int old_count = s->count;
        Entry *old_entries = malloc(old_count * sizeof(Entry));
        if (!old_entries) {
            printf("Ошибка памяти\n");
            return;
        }
        memcpy(old_entries, s->entries, old_count * sizeof(Entry));
        
        // Очищаем состояние (будем заполнять заново)
        s->count = 0;
        
        // Обрабатываем каждую старую запись
        for (int i = 0; i < old_count; i++) {
            int idx = old_entries[i].index;
            double complex amp = old_entries[i].amp;
            int bit = (idx >> qubit) & 1;
            
            // Индексы с битом 0 и 1
            int idx0 = idx;                     // бит qubit = 0
            int idx1 = idx ^ (1 << qubit);      // бит qubit = 1
            
            // Правило преобразования Адамара:
            // |0> -> (|0>+|1>)/√2, |1> -> (|0>-|1>)/√2
            double complex new_amp = amp * factor;
            if (bit == 0) {
                // состояние с битом 0 даёт вклад + в оба
                add_amp(s, idx0,  new_amp);
                add_amp(s, idx1,  new_amp);
            } else {
                // состояние с битом 1 даёт вклад + в idx0, - в idx1
                add_amp(s, idx0,  new_amp);
                add_amp(s, idx1, -new_amp);
            }
        }
        free(old_entries);
    }
}