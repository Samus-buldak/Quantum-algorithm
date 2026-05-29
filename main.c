#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "state.h"
#include "grover.h"

int main() {
    int n, target;
    printf("Введите число кубитов: ");
    scanf("%d", &n);
    int N = 1 << n;
    printf("Введите искомый индекс (0..%d): ", N-1);
    scanf("%d", &target);
    
    freopen("output.txt", "w", stdout);
    srand(time(NULL));
    printf("=== Симуляция квантового алгоритма Гровера ===\n");
    printf("n = %d, N = %d, target = %d\n", n, N, target);
    
    clock_t start = clock();
    SparseState *state = create_state(n);
    if (!state) return 1;
    
    printf("\n1. Начальное состояние |0>:\n");
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));
    
    printf("\n2. Применяем H⊗%d:\n", n);
    apply_hadamard_all(state);
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));
    
    int iterations = (int)(M_PI/4.0 * sqrt(N));
    if (iterations < 1) iterations = 1;
    printf("\n3. Число итераций Гровера: %d\n", iterations);
    
    // Основной цикл Гровера
    for (int i = 0; i < iterations; i++) {
        apply_oracle(state, target);
        apply_diffusion(state);
    }
    
    printf("\n4. Состояние после %d итераций:\n", iterations);
    print_state(state);
    printf("Квадрат нормы: %.10f\n", norm_sq(state));
    
    // Измерение
    double prob;
    int result = measure(state, &prob);
    printf("\n5. Результат измерения:\n");
    printf("   Индекс: %d\n", result);
    printf("   Вероятность: %.6f\n", prob);
    printf("   Ненулевых амплитуд: %d\n", count_nonzero(state));
    
    clock_t end = clock();
    printf("\n6. Время выполнения: %.3f сек\n", (double)(end-start)/CLOCKS_PER_SEC);
    
    // Демонстрация QFT
    printf("\n=== Демонстрация QFT ===\n");
    SparseState *s2 = create_state(n);
    apply_qft(s2, 0);
    printf("QFT от |0⟩:\n");
    print_state(s2);
    apply_qft(s2, 1);
    printf("Обратное QFT (должно быть |0⟩):\n");
    print_state(s2);
    free_state(s2);
    
    free_state(state);
    return 0;
}