#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include "state.h"
#include "grover.h"

int main() {
    int n;
    int target; // искомый индекс

    printf("Введите число кубитов: ");
    scanf("%d", &n);
    int N = 1 << n;
    printf("Введите искомый индекс (от 0 до %d): ", N - 1);
    scanf("%d", &target);
    while (target < 0 || target >= N) {
        printf("Индекс не попадает в диапазон, повторите попытку: ");
        scanf("%d", &target);
    }

    freopen("output.txt", "w", stdout);

    printf("=== Симуляция квантового алгоритма Гровера ===\n");
    printf("Число кубитов: n = %d, размер пространства: N = %d\n", n, N);
    printf("Искомый элемент: %d\n", target);

    // Засекаем время
    clock_t start_time = clock();

    // Создание состояния |0⟩
    SparseState *state = create_state(n);
    if (!state) {
        printf("Ошибка создания состояния\n");
        return 1;
    }

    printf("\n1. Начальное состояние |0>:\n");
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));

    // Применяем H⊗n для создания равномерной суперпозиции
    printf("\n2. Применяем гейт Адамара H⊗%d:\n", n);
    apply_hadamard_all(state);
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));

    // Определяем число итераций Гровера
    int iterations = (int)(M_PI / 4.0 * sqrt(N));
    printf("\n3. Число итераций Гровера: %d\n", iterations);

    // Цикл Гровера: оракул + диффузия
    for (int i = 0; i < iterations; i++) {
        // Оракул: инверсия фазы целевого состояния
        set_amp(state, target, -get_amp(state, target));
        // Оператор диффузии
        apply_diffusion(state);
    }

    printf("\n4. Состояние после %d итераций Гровера:\n", iterations);
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));

    // Измерение
    double probability;
    int result = measure(state, &probability);
    printf("\n5. Результат измерения:\n");
    printf("   Измеренный индекс: %d\n", result);
    printf("   Вероятность: %.6f\n", probability);

    // Подсчёт ненулевых амплитуд
    int nonzero = count_nonzero(state);
    printf("   Ненулевых амплитуд: %d из %d\n", nonzero, N);

    // Время выполнения
    clock_t end_time = clock();
    double cpu_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\n6. Время выполнения симуляции: %.3f секунд\n", cpu_time);

    // Демонстрация QFT и обратного QFT (на отдельном состоянии |0⟩)
    printf("\n=== Дополнительно: демонстрация QFT и обратного QFT ===\n");
    SparseState *s_qft = create_state(n);  // состояние |0⟩
    printf("Начальное состояние |0>:\n");
    print_state(s_qft);

    // Прямое QFT
    apply_qft(s_qft, 0);
    printf("\nПосле прямого QFT:\n");
    print_state(s_qft);

    // Обратное QFT (должно вернуть |0⟩)
    apply_qft(s_qft, 1);
    printf("\nПосле обратного QFT (должно быть исходное |0>):\n");
    print_state(s_qft);

    // Освобождение памяти
    free_state(s_qft);
    free_state(state);
    return 0;
}