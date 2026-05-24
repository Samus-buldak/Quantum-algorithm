#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "state.h"
#include "grover.h"

int main() {
    int n;
    int target; // искомый индекс

    printf("Введите число кубитов: ");
    scanf("%d", &n);
    printf("Введите искомый индекс (от 0 до %.0f): ", pow(2, n));
    scanf("%d", &target);
    while (target < 0 || target >= pow(2, n)){
        printf("индекс не поападает в диапазон, повторите попытку: ");
        scanf("%d", &target);
    }

    freopen("output.txt", "w", stdout);

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

    // Оракул: инверсия фазы целевого состояния
    set_amp(state, target, -get_amp(state, target));   // умножаем на -1

    // Диффузия
    apply_diffusion(state);
    printf("\nПосле оператора диффузии:\n");
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));

    // Проверяем QFT
    SparseState *s = create_state(n);  // |0⟩
    printf("Начальное состояние:\n");
    print_state(s);

    // Применяем QFT для демонстрации
    apply_qft(s, 0);   // прямое QFT
    printf("\nПосле прямого QFT:\n");
    print_state(s);

    // Обратное QFT должно вернуть исходное состояние
    apply_qft(s, 1);   // обратное QFT
    printf("\nПосле обратного QFT (должно быть исходное |0⟩):\n");
    print_state(s);

    free_state(s);
    free_state(state);
    return 0;
}