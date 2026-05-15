"""
Скрипт идентификации упругих характеристик
Использует обертку для вызова исходной программы
Подбирает E1, E2 и G12
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from laminate_wrapper import calculate_laminate_properties


# ВХОДНЫЕ ДАННЫЕ

# Нагрузка, для которой проводим идентификацию
P_LOAD_N = 900 

# НАПРЯЖЕНИЯ

SIGMA_X_FRONT = 11.61399816      # МПа — напряжение вдоль оси x
SIGMA_Y_FRONT = 0.0              # МПа — поперечное напряжение


SIGMA_X_REAR = 11.61388882       # МПа — напряжение вдоль оси x
SIGMA_Y_REAR = 0.0               # МПа — поперечное напряжение

# ЭКСПЕРИМЕНТАЛЬНЫЕ ДЕФОРМАЦИИ 
# Датчик 1 , датчик 2 
EPS_X_EXP_FRONT = 550.0607616 * 10**(-9) 
EPS_Y_EXP_FRONT = 388.8673447 * 10**(-9) 

# Датчик 3 , датчик 4
EPS_X_EXP_REAR = 509.5066568 * 10**(-9) 
EPS_Y_EXP_REAR = 388.9562186 * 10**(-9)  

# ПАРАМЕТРЫ ПАКЕТА 

# Общая толщина пакета
TOTAL_THICKNESS = 8.0  # мм

# Количество слоев
N_LAYERS = 2

# Толщина одного слоя
LAYER_THICKNESS = TOTAL_THICKNESS / N_LAYERS  # мм

# УГЛЫ АРМИРОВАНИЯ
ANGLES_DEG = [30, 30]

# ПАРАМЕТРЫ МАТЕРИАЛА

# Коэффициенты Пуассона
NU12 = 0.3
NU21 = 0.3

# Диапазоны для перебора
E1_RANGE = np.linspace(10000, 90000, 15)   # МПа
E2_RANGE = np.linspace(5000, 50000, 15)    # МПа
G12_RANGE = np.linspace(3000, 5000, 15)    # МПа

# ФУНКЦИЯ РАСПРЕДЕЛЕНИЯ НАПРЯЖЕНИЙ ПО СЛОЯМ

def distribute_stresses(sigma_x, sigma_y, angles_deg, n_layers):
    """
    Распределение глобальных напряжений по слоям
    с учетом углов армирования
    """
    sigma1_list = []
    sigma2_list = []
    tau12_list = []
    
    for angle in angles_deg:
        rad = np.radians(angle)
        c = np.cos(rad)
        s = np.sin(rad)
        c2 = c * c
        s2 = s * s
        cs = c * s
        
        tau_xy = 0.0
        
        sigma1 = sigma_x * c2 + sigma_y * s2 + 2 * tau_xy * cs
        sigma2 = sigma_x * s2 + sigma_y * c2 - 2 * tau_xy * cs
        tau12 = -sigma_x * cs + sigma_y * cs + tau_xy * (c2 - s2)
        
        sigma1_list.append(sigma1)
        sigma2_list.append(sigma2)
        tau12_list.append(tau12)
    
    return sigma1_list, sigma2_list, tau12_list

# ФУНКЦИЯ ОЦЕНКИ ОШИБКИ

def evaluate_error(E1, E2, G12,
                   sigma_x_front, sigma_y_front,
                   sigma_x_rear, sigma_y_rear,
                   eps_x_exp_front, eps_y_exp_front,
                   eps_x_exp_rear, eps_y_exp_rear,
                   angles_deg, n_layers):
    """
    Оценка ошибки для заданных E1, E2, G12
    """
    # Формируем списки свойств
    E1_list = [E1] * n_layers
    E2_list = [E2] * n_layers
    G12_list = [G12] * n_layers
    
    # Распределяем напряжения по слоям для передней поверхности
    sigma1_front, sigma2_front, tau12_front = distribute_stresses(
        sigma_x_front, sigma_y_front, angles_deg, n_layers
    )
    
    # Вызов программы через обертку
    result_front = calculate_laminate_properties(
        E1_list, E2_list, G12_list, NU12, NU21, angles_deg,
        sigma1_front, sigma2_front, tau12_front
    )
    
    sigma1_rear, sigma2_rear, tau12_rear = distribute_stresses(
        sigma_x_rear, sigma_y_rear, angles_deg, n_layers
    )
    
    result_rear = calculate_laminate_properties(
        E1_list, E2_list, G12_list, NU12, NU21, angles_deg,
        sigma1_rear, sigma2_rear, tau12_rear
    )
    
    if result_front['eps_x'] is None or result_rear['eps_x'] is None:
        return float('inf'), None, None
    
    # Ошибка: сумма квадратов отклонений
    # Взвешиваем ошибки, чтобы учесть разный масштаб деформаций
    error = (
        ((result_front['eps_x'] - eps_x_exp_front) / eps_x_exp_front)**2 +
        ((result_front['eps_y'] - eps_y_exp_front) / eps_y_exp_front)**2 +
        ((result_rear['eps_x'] - eps_x_exp_rear) / eps_x_exp_rear)**2 +
        ((result_rear['eps_y'] - eps_y_exp_rear) / eps_y_exp_rear)**2
    )
    
    return error, result_front, result_rear

# ФУНКЦИЯ ДЛЯ ПОСТРОЕНИЯ СЕЧЕНИЙ КАРТЫ ОШИБОК

def plot_error_slices(best_E1, best_E2, best_G12, best_error,
                      E1_range, E2_range, G12_range,
                      sigma_x_front, sigma_y_front,
                      sigma_x_rear, sigma_y_rear,
                      eps_x_exp_front, eps_y_exp_front,
                      eps_x_exp_rear, eps_y_exp_rear,
                      angles_deg, n_layers):
    """
    Построение сечений карты ошибок для визуализации влияния каждого параметра
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    print("\nПостроение сечений карты ошибок...")
    
    error_grid_E1E2 = np.zeros((len(E1_range), len(E2_range)))
    for i, E1 in enumerate(E1_range):
        for j, E2 in enumerate(E2_range):
            error, _, _ = evaluate_error(
                E1, E2, best_G12,
                sigma_x_front, sigma_y_front,
                sigma_x_rear, sigma_y_rear,
                eps_x_exp_front, eps_y_exp_front,
                eps_x_exp_rear, eps_y_exp_rear,
                angles_deg, n_layers
            )
            error_grid_E1E2[i, j] = error
    
    X, Y = np.meshgrid(E2_range, E1_range)
    log_error = np.log10(error_grid_E1E2 + 1e-10)
    
    contour1 = axes[0].contourf(X, Y, log_error, levels=20, cmap='viridis')
    axes[0].scatter(best_E2, best_E1, color='red', marker='*', s=200,
                    label=f'Оптимум')
    axes[0].set_xlabel('E2, МПа')
    axes[0].set_ylabel('E1, МПа')
    axes[0].set_title(f'Сечение при G12 = {best_G12:.0f} МПа')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    plt.colorbar(contour1, ax=axes[0], label='log10(error)')
    
    error_grid_E1G12 = np.zeros((len(E1_range), len(G12_range)))
    for i, E1 in enumerate(E1_range):
        for k, G12 in enumerate(G12_range):
            error, _, _ = evaluate_error(
                E1, best_E2, G12,
                sigma_x_front, sigma_y_front,
                sigma_x_rear, sigma_y_rear,
                eps_x_exp_front, eps_y_exp_front,
                eps_x_exp_rear, eps_y_exp_rear,
                angles_deg, n_layers
            )
            error_grid_E1G12[i, k] = error
    
    X2, Y2 = np.meshgrid(G12_range, E1_range)
    log_error2 = np.log10(error_grid_E1G12 + 1e-10)
    
    contour2 = axes[1].contourf(X2, Y2, log_error2, levels=20, cmap='viridis')
    axes[1].scatter(best_G12, best_E1, color='red', marker='*', s=200,
                    label=f'Оптимум')
    axes[1].set_xlabel('G12, МПа')
    axes[1].set_ylabel('E1, МПа')
    axes[1].set_title(f'Сечение при E2 = {best_E2:.0f} МПа')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    plt.colorbar(contour2, ax=axes[1], label='log10(error)')
    
    error_grid_E2G12 = np.zeros((len(E2_range), len(G12_range)))
    for j, E2 in enumerate(E2_range):
        for k, G12 in enumerate(G12_range):
            error, _, _ = evaluate_error(
                best_E1, E2, G12,
                sigma_x_front, sigma_y_front,
                sigma_x_rear, sigma_y_rear,
                eps_x_exp_front, eps_y_exp_front,
                eps_x_exp_rear, eps_y_exp_rear,
                angles_deg, n_layers
            )
            error_grid_E2G12[j, k] = error
    
    X3, Y3 = np.meshgrid(G12_range, E2_range)
    log_error3 = np.log10(error_grid_E2G12 + 1e-10)
    
    contour3 = axes[2].contourf(X3, Y3, log_error3, levels=20, cmap='viridis')
    axes[2].scatter(best_G12, best_E2, color='red', marker='*', s=200,
                    label=f'Оптимум')
    axes[2].set_xlabel('G12, МПа')
    axes[2].set_ylabel('E2, МПа')
    axes[2].set_title(f'Сечение при E1 = {best_E1:.0f} МПа')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    plt.colorbar(contour3, ax=axes[2], label='log10(error)')
    
    plt.tight_layout()
    plt.savefig('identification_slices.png', dpi=150)
    plt.show()


def identify_properties():
    """
    Идентификация свойств методом перебора по трем параметрам: E1, E2, G12
    """
    print("=" * 80)
    print("ИДЕНТИФИКАЦИЯ УПРУГИХ ХАРАКТЕРИСТИК КОМПОЗИТНОГО ПАКЕТА")
    print("=" * 80)
    print()
    print("Исходные данные:")
    print(f"  Нагрузка: {P_LOAD_N} Н ({P_LOAD_N/9.81:.1f} кг)")
    print(f"  Напряжение на передней поверхности: σx = {SIGMA_X_FRONT:.4f} МПа")
    print(f"  Напряжение на задней поверхности: σx = {SIGMA_X_REAR:.4f} МПа")
    print(f"  Экспериментальные деформации (передняя):")
    print(f"    εx = {EPS_X_EXP_FRONT:.2e}, εy = {EPS_Y_EXP_FRONT:.2e}")
    print(f"  Экспериментальные деформации (задняя):")
    print(f"    εx = {EPS_X_EXP_REAR:.2e}, εy = {EPS_Y_EXP_REAR:.2e}")
    print(f"  Углы армирования: {ANGLES_DEG}")
    print(f"  Количество слоев: {N_LAYERS}")
    print(f"  Толщина пакета: {TOTAL_THICKNESS} мм")
    print()
    
    print("Диапазоны перебора:")
    print(f"  E1: {E1_RANGE[0]:.0f} ... {E1_RANGE[-1]:.0f} МПа ({len(E1_RANGE)} значений)")
    print(f"  E2: {E2_RANGE[0]:.0f} ... {E2_RANGE[-1]:.0f} МПа ({len(E2_RANGE)} значений)")
    print(f"  G12: {G12_RANGE[0]:.0f} ... {G12_RANGE[-1]:.0f} МПа ({len(G12_RANGE)} значений)")
    
    total_combinations = len(E1_RANGE) * len(E2_RANGE) * len(G12_RANGE)
    print(f"\nВсего комбинаций: {total_combinations}")
    print()
    
    # Поиск оптимальных значений
    best_error = float('inf')
    best_E1 = None
    best_E2 = None
    best_G12 = None
    best_results = None
    
    # Счетчик для прогресса
    counter = 0
    
    print("Выполняется перебор параметров...")
    print()
    
    # Вложенные циклы перебора
    for i, E1 in enumerate(E1_RANGE):
        for j, E2 in enumerate(E2_RANGE):
            for k, G12 in enumerate(G12_RANGE):
                counter += 1
                
                # Выводим прогресс каждые 100 итераций
                if counter % 100 == 0 or counter == total_combinations:
                    print(f"  Прогресс: {counter}/{total_combinations} ({100*counter/total_combinations:.1f}%)")
                
                error, res_front, res_rear = evaluate_error(
                    E1, E2, G12,
                    SIGMA_X_FRONT, SIGMA_Y_FRONT,
                    SIGMA_X_REAR, SIGMA_Y_REAR,
                    EPS_X_EXP_FRONT, EPS_Y_EXP_FRONT,
                    EPS_X_EXP_REAR, EPS_Y_EXP_REAR,
                    ANGLES_DEG, N_LAYERS
                )
                
                if error < best_error:
                    best_error = error
                    best_E1 = E1
                    best_E2 = E2
                    best_G12 = G12
                    best_results = (res_front, res_rear)
                    
                    print(f"\n  ★ НОВЫЙ ЛУЧШИЙ РЕЗУЛЬТАТ ★")
                    print(f"    E1 = {E1:.0f} МПа, E2 = {E2:.0f} МПа, G12 = {G12:.0f} МПа")
                    print(f"    Ошибка = {error:.4e}")
                    print()
    
    print()
    print("=" * 80)
    print("РЕЗУЛЬТАТЫ ИДЕНТИФИКАЦИИ")
    print("=" * 80)
    print()
    print("Оптимальные упругие характеристики:")
    print(f"  E1 = {best_E1:.0f} МПа ({best_E1/1000:.1f} ГПа)")
    print(f"  E2 = {best_E2:.0f} МПа ({best_E2/1000:.1f} ГПа)")
    print(f"  G12 = {best_G12:.0f} МПа ({best_G12/1000:.1f} ГПа)")
    print(f"  Коэффициент Пуассона ν12 = {NU12}")
    print(f"  Коэффициент Пуассона ν21 = {NU21}")
    print(f"  Минимальная ошибка: {best_error:.4e}")
    print()
    
    if best_results:
        res_front, res_rear = best_results
        print("Сравнение деформаций:")
        print()
        print("  Передняя поверхность:")
        print(f"    Расчетные: εx = {res_front['eps_x']:.6e}, εy = {res_front['eps_y']:.6e}")
        print(f"    Эксперимент: εx = {EPS_X_EXP_FRONT:.6e}, εy = {EPS_Y_EXP_FRONT:.6e}")
        print(f"    Относительное отклонение:")
        print(f"      Δεx = {abs(res_front['eps_x'] - EPS_X_EXP_FRONT) / abs(EPS_X_EXP_FRONT) * 100:.2f}%")
        print(f"      Δεy = {abs(res_front['eps_y'] - EPS_Y_EXP_FRONT) / abs(EPS_Y_EXP_FRONT) * 100:.2f}%")
        print()
        print("  Задняя поверхность:")
        print(f"    Расчетные: εx = {res_rear['eps_x']:.6e}, εy = {res_rear['eps_y']:.6e}")
        print(f"    Эксперимент: εx = {EPS_X_EXP_REAR:.6e}, εy = {EPS_Y_EXP_REAR:.6e}")
        print(f"    Относительное отклонение:")
        print(f"      Δεx = {abs(res_rear['eps_x'] - EPS_X_EXP_REAR) / abs(EPS_X_EXP_REAR) * 100:.2f}%")
        print(f"      Δεy = {abs(res_rear['eps_y'] - EPS_Y_EXP_REAR) / abs(EPS_Y_EXP_REAR) * 100:.2f}%")
    
    # Построение сечений карты ошибок
    print()
    print("Построение графиков...")
    plot_error_slices(
        best_E1, best_E2, best_G12, best_error,
        E1_RANGE, E2_RANGE, G12_RANGE,
        SIGMA_X_FRONT, SIGMA_Y_FRONT,
        SIGMA_X_REAR, SIGMA_Y_REAR,
        EPS_X_EXP_FRONT, EPS_Y_EXP_FRONT,
        EPS_X_EXP_REAR, EPS_Y_EXP_REAR,
        ANGLES_DEG, N_LAYERS
    )
    
    return {
        'E1': best_E1,
        'E2': best_E2,
        'G12': best_G12,
        'nu12': NU12,
        'nu21': NU21,
        'error': best_error,
        'results_front': best_results[0] if best_results else None,
        'results_rear': best_results[1] if best_results else None
    }


if __name__ == "__main__":
    results = identify_properties()
    
    # Вывод итоговых значений для копирования в predictor.py
    print()
    print("=" * 80)
    print("ЗНАЧЕНИЯ ДЛЯ ПОДСТАНОВКИ В predictor.py")
    print("=" * 80)
    print(f"""
# Идентифицированные свойства (из identification.py)
E1_IDENTIFIED = {results['E1']:.0f}   # МПа
E2_IDENTIFIED = {results['E2']:.0f}   # МПа
G12_IDENTIFIED = {results['G12']:.0f} # МПа
NU12 = {results['nu12']}
NU21 = {results['nu21']}
    """)