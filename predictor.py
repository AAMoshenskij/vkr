"""
Модуль прогнозирования поведения протеза стопы
Использует идентифицированные свойства для прогноза при различных нагрузках
"""
import numpy as np
import matplotlib.pyplot as plt
from laminate_wrapper import calculate_laminate_properties


class ProsthesisPredictor:
    """
    Система прогнозирования поведения протеза стопы
    """
    
    def __init__(self, E1, E2, G12, nu12, nu21, angles_deg, n_layers, 
                 strength_tension=1500, strength_compression=1200, 
                 strength_shear=80, strength_transverse=50):
        """
        Инициализация с идентифицированными свойствами
        
        Параметры:
        -----------
        E1, E2, G12 : float
            Упругие характеристики материала (МПа)
        nu12, nu21 : float
            Коэффициенты Пуассона
        angles_deg : list
            Углы армирования слоев
        n_layers : int
            Количество слоев
        strength_tension, strength_compression : float
            Пределы прочности вдоль волокон (МПа)
        strength_shear : float
            Предел прочности при сдвиге (МПа)
        strength_transverse : float
            Предел прочности поперек волокон (МПа)
        """
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        self.nu12 = nu12
        self.nu21 = nu21
        self.angles_deg = angles_deg
        self.n_layers = n_layers
        
        # Пределы прочности (МПа)
        self.strength_tension = strength_tension      # растяжение вдоль
        self.strength_compression = strength_compression  # сжатие вдоль
        self.strength_shear = strength_shear          # сдвиг
        self.strength_transverse = strength_transverse    # поперек волокон
        
    def _distribute_stresses(self, sigma_x, sigma_y, tau_xy=0.0):
        """
        Распределение глобальных напряжений по слоям
        с полным преобразованием тензора напряжений
        
        Вход:
            sigma_x, sigma_y, tau_xy — глобальные напряжения (МПа)
        
        Выход:
            sigma1_list, sigma2_list, tau12_list — напряжения в каждом слое
        """
        sigma1_list = []
        sigma2_list = []
        tau12_list = []
        
        for angle in self.angles_deg:
            rad = np.radians(angle)
            c = np.cos(rad)
            s = np.sin(rad)
            c2 = c * c
            s2 = s * s
            cs = c * s
            
            # Преобразование напряжений в локальную систему координат слоя
            sigma1 = sigma_x * c2 + sigma_y * s2 + 2 * tau_xy * cs
            sigma2 = sigma_x * s2 + sigma_y * c2 - 2 * tau_xy * cs
            tau12 = -sigma_x * cs + sigma_y * cs + tau_xy * (c2 - s2)
            
            sigma1_list.append(sigma1)
            sigma2_list.append(sigma2)
            tau12_list.append(tau12)
        
        return sigma1_list, sigma2_list, tau12_list
    
    def predict_for_load(self, sigma_x_front, sigma_y_front, tau_xy_front,
                           sigma_x_rear, sigma_y_rear, tau_xy_rear):
        """
        Прогнозирование деформаций для заданных напряжений
        """
        # Свойства для всех слоев одинаковые
        E1_list = [self.E1] * self.n_layers
        E2_list = [self.E2] * self.n_layers
        G12_list = [self.G12] * self.n_layers
        
        # Расчет для передней поверхности
        sigma1_front, sigma2_front, tau12_front = self._distribute_stresses(
            sigma_x_front, sigma_y_front, tau_xy_front
        )
        result_front = calculate_laminate_properties(
            E1_list, E2_list, G12_list, self.nu12, self.nu21, self.angles_deg,
            sigma1_front, sigma2_front, tau12_front
        )
        
        # Расчет для задней поверхности
        sigma1_rear, sigma2_rear, tau12_rear = self._distribute_stresses(
            sigma_x_rear, sigma_y_rear, tau_xy_rear
        )
        result_rear = calculate_laminate_properties(
            E1_list, E2_list, G12_list, self.nu12, self.nu21, self.angles_deg,
            sigma1_rear, sigma2_rear, tau12_rear
        )
        
        return {
            'front': result_front,
            'rear': result_rear,
            'layer_stresses_front': (sigma1_front, sigma2_front, tau12_front),
            'layer_stresses_rear': (sigma1_rear, sigma2_rear, tau12_rear)
        }
    
    def check_safety_tsai_wu(self, sigma1, sigma2, tau12):
        """
        Проверка прочности по критерию Цая-Ву
        
        Критерий: F1·σ1 + F2·σ2 + F11·σ1² + F22·σ2² + F66·τ12² + 2F12·σ1·σ2 ≤ 1
        
        Параметры:
        -----------
        sigma1, sigma2, tau12 : float
            Напряжения в слое (МПа)
            
        Возвращает:
        -----------
        float : индекс Цая-Ву (<1 — безопасно, >1 — разрушение)
        """
        # Коэффициенты прочности (для углепластика)
        Xt = self.strength_tension      # прочность при растяжении вдоль
        Xc = self.strength_compression  # прочность при сжатии вдоль
        Yt = self.strength_transverse   # прочность при растяжении поперек
        Yc = self.strength_transverse   # прочность при сжатии поперек
        S = self.strength_shear         # прочность при сдвиге
        
        # Параметры Tsai-Wu
        F1 = 1/Xt - 1/Xc
        F2 = 1/Yt - 1/Yc
        F11 = 1/(Xt * Xc)
        F22 = 1/(Yt * Yc)
        F66 = 1/(S * S)
        
        # Взаимодействующий член
        F12 = -0.5 * np.sqrt(F11 * F22)
        
        # Вычисление индекса
        tsai_wu_index = (F1 * sigma1 + F2 * sigma2 + 
                         F11 * sigma1**2 + F22 * sigma2**2 + 
                         F66 * tau12**2 + 2 * F12 * sigma1 * sigma2)
        
        return tsai_wu_index
    
    def check_safety_max_stress(self, sigma1, sigma2, tau12):
        """
        Проверка прочности по критерию максимальных напряжений
        
        Возвращает:
        -----------
        dict : результаты проверки для каждого типа напряжений
        """
        # Коэффициенты запаса по каждому типу напряжений
        if sigma1 > 0:
            sf_1 = self.strength_tension / sigma1 if sigma1 > 0 else float('inf')
        else:
            sf_1 = abs(self.strength_compression / sigma1) if sigma1 < 0 else float('inf')
        
        if sigma2 > 0:
            sf_2 = self.strength_transverse / sigma2 if sigma2 > 0 else float('inf')
        else:
            sf_2 = abs(self.strength_transverse / sigma2) if sigma2 < 0 else float('inf')
        
        sf_12 = self.strength_shear / abs(tau12) if abs(tau12) > 0 else float('inf')
        
        # Минимальный коэффициент запаса
        min_sf = min(sf_1, sf_2, sf_12)
        
        return {
            'sf_1': sf_1,
            'sf_2': sf_2,
            'sf_12': sf_12,
            'min_sf': min_sf
        }
    
    def analyze_layer_safety(self, sigma1_list, sigma2_list, tau12_list, layer_names=None):
        """
        Анализ безопасности каждого слоя
        """
        results = []
        
        for i, (s1, s2, t12) in enumerate(zip(sigma1_list, sigma2_list, tau12_list)):
            # Критерий максимальных напряжений
            max_stress_result = self.check_safety_max_stress(s1, s2, t12)
            
            # Критерий Цая-Ву
            tsai_wu = self.check_safety_tsai_wu(s1, s2, t12)
            
            # Статус слоя
            if tsai_wu < 0.7:
                status = "Безопасно"
                color = "green"
            elif tsai_wu < 0.9:
                status = "Допустимо"
                color = "yellow"
            elif tsai_wu < 1.0:
                status = "Критично"
                color = "orange"
            else:
                status = "Разрушение"
                color = "red"
            
            results.append({
                'layer': i,
                'angle': self.angles_deg[i],
                'sigma1': s1,
                'sigma2': s2,
                'tau12': t12,
                'max_stress_sf': max_stress_result['min_sf'],
                'tsai_wu_index': tsai_wu,
                'status': status,
                'color': color
            })
        
        return results
    
    def predict_with_layer_analysis(self, sigma_x_front, sigma_y_front, tau_xy_front,
                                sigma_x_rear, sigma_y_rear, tau_xy_rear):
        """
        Прогнозирование с полным анализом по слоям
        """
        # Получаем деформации и напряжения в слоях
        results = self.predict_for_load(
            sigma_x_front, sigma_y_front, tau_xy_front,
            sigma_x_rear, sigma_y_rear, tau_xy_rear
        )
        
        # Анализ безопасности слоев на передней поверхности
        s1_front, s2_front, t12_front = results['layer_stresses_front']
        layer_analysis_front = self.analyze_layer_safety(s1_front, s2_front, t12_front)
        
        # Анализ безопасности слоев на задней поверхности
        s1_rear, s2_rear, t12_rear = results['layer_stresses_rear']
        layer_analysis_rear = self.analyze_layer_safety(s1_rear, s2_rear, t12_rear)
        
        # Общая оценка конструкции (по наиболее опасному слою)
        all_tsai_wu = [l['tsai_wu_index'] for l in layer_analysis_front + layer_analysis_rear]
        max_tsai_wu = max(all_tsai_wu)
        worst_layer_index = np.argmax(all_tsai_wu)
        worst_layer = (layer_analysis_front + layer_analysis_rear)[worst_layer_index]
        
        if max_tsai_wu < 0.7:
            overall_status = "КОНСТРУКЦИЯ БЕЗОПАСНА"
            overall_recommendation = "Все слои работают в упругой области. Конструкция пригодна для эксплуатации."
        elif max_tsai_wu < 0.9:
            overall_status = "ДОПУСТИМО, НО ТРЕБУЕТ КОНТРОЛЯ"
            overall_recommendation = "Некоторые слои работают вблизи предела. Рекомендуется снизить нагрузку или усилить конструкцию."
        elif max_tsai_wu < 1.0:
            overall_status = "КРИТИЧЕСКОЕ СОСТОЯНИЕ"
            overall_recommendation = "Конструкция работает на пределе. Эксплуатация возможна только при сниженных нагрузках."
        else:
            overall_status = "ПРОГНОЗИРУЕТСЯ РАЗРУШЕНИЕ"
            overall_recommendation = "Критерий Цая-Ву превышен. Конструкция не пригодна для эксплуатации."
        
        return {
            'deformations': results,
            'layer_analysis_front': layer_analysis_front,
            'layer_analysis_rear': layer_analysis_rear,
            'max_tsai_wu': max_tsai_wu,
            'worst_layer': worst_layer,
            'overall_status': overall_status,
            'overall_recommendation': overall_recommendation
        }

    
    def generate_detailed_report(self, load_kg, sigma_x_front, sigma_x_rear,
                              sigma_y_front=0, sigma_y_rear=0,
                              tau_xy_front=0, tau_xy_rear=0):
        """
        Отчет с анализом по слоям
        """
        print("=" * 80)
        print("СИСТЕМА ПРОГНОЗИРОВАНИЯ ПОВЕДЕНИЯ ПРОТЕЗА СТОПЫ")
        print("=" * 80)
        
        print("\n" + "─" * 80)
        print("ВХОДНЫЕ ДАННЫЕ")
        print("─" * 80)
        print(f"  Нагрузка: {load_kg:.1f} кг ({load_kg * 9.81:.0f} Н)")
        print(f"  Напряжения на передней поверхности: σx = {sigma_x_front:.2f} МПа, σy = {sigma_y_front:.2f} МПа")
        print(f"  Напряжения на задней поверхности:  σx = {sigma_x_rear:.2f} МПа, σy = {sigma_y_rear:.2f} МПа")
        
        print("\n" + "─" * 80)
        print("ИДЕНТИФИЦИРОВАННЫЕ СВОЙСТВА МАТЕРИАЛА")
        print("─" * 80)
        print(f"  E1 = {self.E1:.0f} МПа ({self.E1/1000:.1f} ГПа)")
        print(f"  E2 = {self.E2:.0f} МПа ({self.E2/1000:.1f} ГПа)")
        print(f"  G12 = {self.G12:.0f} МПа ({self.G12/1000:.1f} ГПа)")
        print(f"  ν12 = {self.nu12:.3f}, ν21 = {self.nu21:.3f}")
        print(f"  Структура: {self.n_layers} слоев, углы {self.angles_deg}")
        
        print("\n" + "─" * 80)
        print("ПРОЧНОСТНЫЕ ХАРАКТЕРИСТИКИ МАТЕРИАЛА")
        print("─" * 80)
        print(f"  Предел прочности при растяжении вдоль: {self.strength_tension:.0f} МПа")
        print(f"  Предел прочности при сжатии вдоль:    {self.strength_compression:.0f} МПа")
        print(f"  Предел прочности поперек волокон:     {self.strength_transverse:.0f} МПа")
        print(f"  Предел прочности при сдвиге:          {self.strength_shear:.0f} МПа")
        
        # Прогноз с анализом по слоям
        analysis = self.predict_with_layer_analysis(
            sigma_x_front, sigma_y_front, tau_xy_front,
            sigma_x_rear, sigma_y_rear, tau_xy_rear
        )
        
        # Деформации
        print("\n" + "─" * 80)
        print("ПРОГНОЗИРУЕМЫЕ ДЕФОРМАЦИИ")
        print("─" * 80)
        res_front = analysis['deformations']['front']
        res_rear = analysis['deformations']['rear']
        print(f"  Передняя поверхность: εx = {res_front['eps_x']:.6e} ({res_front['eps_x']*1e6:.1f} мкм/м)")
        print(f"                        εy = {res_front['eps_y']:.6e} ({res_front['eps_y']*1e6:.1f} мкм/м)")
        print(f"  Задняя поверхность:   εx = {res_rear['eps_x']:.6e} ({res_rear['eps_x']*1e6:.1f} мкм/м)")
        print(f"                        εy = {res_rear['eps_y']:.6e} ({res_rear['eps_y']*1e6:.1f} мкм/м)")
        
        # Анализ слоев (передняя поверхность)
        print("\n" + "─" * 80)
        print("АНАЛИЗ СЛОЕВ (ПЕРЕДНЯЯ ПОВЕРХНОСТЬ)")
        print("─" * 80)
        print(f"{'Слой':>4} | {'Угол':>5} | {'σ1, МПа':>10} | {'σ2, МПа':>10} | {'τ12, МПа':>10} | {'Крит. Ц-В':>10} | {'Статус':>12}")
        print("-" * 80)
        for layer in analysis['layer_analysis_front']:
            print(f"{layer['layer']:>4} | {layer['angle']:>5}° | {layer['sigma1']:>10.2f} | {layer['sigma2']:>10.2f} | {layer['tau12']:>10.2f} | {layer['tsai_wu_index']:>10.3f} | {layer['status']:>12}")
        
        # Анализ слоев (задняя поверхность)
        print("\n" + "─" * 80)
        print("АНАЛИЗ СЛОЕВ (ЗАДНЯЯ ПОВЕРХНОСТЬ)")
        print("─" * 80)
        print(f"{'Слой':>4} | {'Угол':>5} | {'σ1, МПа':>10} | {'σ2, МПа':>10} | {'τ12, МПа':>10} | {'Крит. Ц-В':>10} | {'Статус':>12}")
        print("-" * 80)
        for layer in analysis['layer_analysis_rear']:
            print(f"{layer['layer']:>4} | {layer['angle']:>5}° | {layer['sigma1']:>10.2f} | {layer['sigma2']:>10.2f} | {layer['tau12']:>10.2f} | {layer['tsai_wu_index']:>10.3f} | {layer['status']:>12}")
        
        # Общий вердикт
        print("\n" + "─" * 80)
        print("ОБЩИЙ ВЕРДИКТ")
        print("─" * 80)
        print(f"  Максимальный индекс Цая-Ву: {analysis['max_tsai_wu']:.3f}")
        print(f"  Наиболее опасный слой: слой {analysis['worst_layer']['layer']} (угол {analysis['worst_layer']['angle']}°)")
        print(f"    Напряжения в слое: σ1 = {analysis['worst_layer']['sigma1']:.2f} МПа, "
            f"σ2 = {analysis['worst_layer']['sigma2']:.2f} МПа, τ12 = {analysis['worst_layer']['tau12']:.2f} МПа")
        print(f"    Коэффициент запаса (по макс. напряжениям): {analysis['worst_layer']['max_stress_sf']:.2f}")
        
        print(f"\n  {analysis['overall_status']}")
        print(f"  {analysis['overall_recommendation']}")
        
        print("\n" + "=" * 80)
        
        return analysis
    
    def plot_layer_stresses(self, sigma_x_front, sigma_x_rear, save_path=None):
        """
        Визуализация распределения напряжений по слоям
        """
        # Получаем напряжения в слоях
        s1_front, s2_front, t12_front = self._distribute_stresses(sigma_x_front, 0, 0)
        s1_rear, s2_rear, t12_rear = self._distribute_stresses(sigma_x_rear, 0, 0)
        
        layers = range(self.n_layers)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # σ1 по слоям
        axes[0, 0].bar(layers, s1_front, alpha=0.7, label='Передняя поверхность', color='blue')
        axes[0, 0].bar(layers, s1_rear, alpha=0.7, label='Задняя поверхность', color='red')
        axes[0, 0].axhline(y=self.strength_tension, color='g', linestyle='--', label='Предел растяжения')
        axes[0, 0].axhline(y=-self.strength_compression, color='orange', linestyle='--', label='Предел сжатия')
        axes[0, 0].set_xlabel('Номер слоя')
        axes[0, 0].set_ylabel('σ1, МПа')
        axes[0, 0].set_title('Напряжения вдоль волокон')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # σ2 по слоям
        axes[0, 1].bar(layers, s2_front, alpha=0.7, label='Передняя поверхность', color='blue')
        axes[0, 1].bar(layers, s2_rear, alpha=0.7, label='Задняя поверхность', color='red')
        axes[0, 1].axhline(y=self.strength_transverse, color='g', linestyle='--', label='Предел поперечный')
        axes[0, 1].axhline(y=-self.strength_transverse, color='orange', linestyle='--')
        axes[0, 1].set_xlabel('Номер слоя')
        axes[0, 1].set_ylabel('σ2, МПа')
        axes[0, 1].set_title('Напряжения поперек волокон')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # τ12 по слоям
        axes[1, 0].bar(layers, t12_front, alpha=0.7, label='Передняя поверхность', color='blue')
        axes[1, 0].bar(layers, t12_rear, alpha=0.7, label='Задняя поверхность', color='red')
        axes[1, 0].axhline(y=self.strength_shear, color='g', linestyle='--', label='Предел сдвига')
        axes[1, 0].axhline(y=-self.strength_shear, color='orange', linestyle='--')
        axes[1, 0].set_xlabel('Номер слоя')
        axes[1, 0].set_ylabel('τ12, МПа')
        axes[1, 0].set_title('Касательные напряжения')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Индекс Цая-Ву
        tsai_wu_front = [self.check_safety_tsai_wu(s1, s2, t12) 
                         for s1, s2, t12 in zip(s1_front, s2_front, t12_front)]
        tsai_wu_rear = [self.check_safety_tsai_wu(s1, s2, t12) 
                        for s1, s2, t12 in zip(s1_rear, s2_rear, t12_rear)]
        
        axes[1, 1].bar(layers, tsai_wu_front, alpha=0.7, label='Передняя поверхность', color='blue')
        axes[1, 1].bar(layers, tsai_wu_rear, alpha=0.7, label='Задняя поверхность', color='red')
        axes[1, 1].axhline(y=1.0, color='red', linestyle='--', label='Критический уровень')
        axes[1, 1].axhline(y=0.7, color='orange', linestyle='--', label='Зона внимания')
        axes[1, 1].set_xlabel('Номер слоя')
        axes[1, 1].set_ylabel('Индекс Цая-Ву')
        axes[1, 1].set_title('Критерий Цая-Ву (≤1 — безопасно)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150)
        
        plt.show()


def main():
    """
    Система прогнозирования
    """
    # Идентифицированные свойства
    E1_IDENTIFIED = 61429  # МПа
    E2_IDENTIFIED = 40357  # МПа
    G12_IDENTIFIED = 3000  # МПа
    NU12 = 0.3
    NU21 = 0.3
    
    # Параметры пакета
    N_LAYERS = 2
    ANGLES_DEG = [90, 0]
    
    # Прочностные характеристики (МПа)
    STRENGTH_TENSION = 1500      # растяжение вдоль волокон
    STRENGTH_COMPRESSION = 1200  # сжатие вдоль волокон
    STRENGTH_SHEAR = 80          # сдвиг
    STRENGTH_TRANSVERSE = 50     # поперек волокон
    
    # Расчет напряжений по формуле прямой балки
    L = 4.0      # мм — плечо
    y = 4.0      # мм — расстояние от нейтральной оси до поверхности
    I = 1621.33  # мм⁴ — момент инерции сечения
    
    # Создание системы прогнозирования
    predictor = ProsthesisPredictor(
        E1_IDENTIFIED, E2_IDENTIFIED, G12_IDENTIFIED, 
        NU12, NU21, ANGLES_DEG, N_LAYERS,
        STRENGTH_TENSION, STRENGTH_COMPRESSION, 
        STRENGTH_SHEAR, STRENGTH_TRANSVERSE
    )
    
    # Функция пересчета нагрузки в напряжение
    def load_to_stress(load_kg):
        """
        Пересчет нагрузки в напряжение по формуле прямой балки
        sigma = F * L * y / I
        """
        F = load_kg * 9.81  # перевод в Ньютоны
        sigma = F * L * y / I  # МПа
        return sigma
    
    # Прогноз для экспериментальных нагрузок

    
    # Данные для прогноза (из эксперимента)
    loads_kg_exp = [30, 60, 90, 120, 135]
    sigma_x_fronts_exp = [3.872, 7.742, 11.614, 15.487, 17.420]
    sigma_x_rears_exp = [3.871, 7.742, 11.614, 15.485, 17.422]
    
    print("\n" + "=" * 80)
    print("ЧАСТЬ 1: ПРОГНОЗИРОВАНИЕ ДЛЯ ЭКСПЕРИМЕНТАЛЬНЫХ НАГРУЗОК")
    print("=" * 80)
    
    for i, load in enumerate(loads_kg_exp):
        print(f"\n{'='*80}")
        print(f"НАГРУЗКА: {load} кг (экспериментальная)")
        print(f"{'='*80}")
        
        predictor.generate_detailed_report(
            load,
            sigma_x_fronts_exp[i], sigma_x_rears_exp[i],
            sigma_y_front=0, sigma_y_rear=0
        )
    
    # Прогнозирование для нагрузок за пределами эксперимента
    
    print("\n" + "=" * 80)
    print("ЧАСТЬ 2: РАСШИРЕННОЕ ПРОГНОЗИРОВАНИЕ (ЗА ПРЕДЕЛАМИ ЭКСПЕРИМЕНТА)")
    print("=" * 80)
    print("\nВ данном разделе выполняется прогнозирование поведения конструкции")
    print("при нагрузках, превышающих экспериментальный диапазон. Расчет")
    print("напряжений выполняется по формуле прямой балки: σ = F·L·y / I")
    print(f"Параметры: L = {L} мм, y = {y} мм, I = {I:.2f} мм⁴")
    
    # Диапазон нагрузок для расширенного прогнозирования
    loads_kg_ext = [150, 180, 210, 240, 270, 300, 350, 400, 450, 500]
    
    # Сбор результатов для анализа
    extended_results = []
    
    for load in loads_kg_ext:
        # Расчет напряжений по формуле прямой балки
        sigma = load_to_stress(load)
        
        print(f"\n{'='*80}")
        print(f"НАГРУЗКА: {load} кг (прогноз)")
        print(f"Расчетное напряжение: σ = {sigma:.2f} МПа")
        print(f"{'='*80}")
        
        # Выполняем прогнозирование
        analysis = predictor.generate_detailed_report(
            load,
            sigma, sigma,
            sigma_y_front=0, sigma_y_rear=0
        )
        
        # Сохраняем результаты для сводной таблицы
        extended_results.append({
            'load_kg': load,
            'sigma': sigma,
            'max_tsai_wu': analysis['max_tsai_wu'],
            'status': analysis['overall_status'],
            'worst_layer_angle': analysis['worst_layer']['angle'],
            'worst_layer_stress': analysis['worst_layer']['sigma2']
        })
    
    # Определение предельной нагрузки
    
    print("\n" + "=" * 80)
    print("ЧАСТЬ 3: СВОДНЫЙ АНАЛИЗ И ОПРЕДЕЛЕНИЕ ПРЕДЕЛЬНОЙ НАГРУЗКИ")
    print("=" * 80)
    
    print("\n" + "─" * 80)
    print("СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ РАСШИРЕННОГО ПРОГНОЗИРОВАНИЯ")
    print("─" * 80)
    print(f"{'Нагрузка, кг':>12} | {'Напряжение, МПа':>15} | {'Индекс Цая-Ву':>14} | {'Статус':>25}")
    print("-" * 80)
    
    for res in extended_results:
        print(f"{res['load_kg']:>12} | {res['sigma']:>15.2f} | {res['max_tsai_wu']:>14.3f} | {res['status']:>25}")
    
    # Определение нагрузок, при которых достигаются различные состояния
    print("\n" + "─" * 80)
    print("АНАЛИЗ ПРЕДЕЛЬНЫХ СОСТОЯНИЙ")
    print("─" * 80)
    
    # Поиск нагрузок, соответствующих граничным значениям индекса Цая-Ву
    
    load_at_07 = None
    load_at_09 = None
    load_at_10 = None
    
    for i in range(len(extended_results) - 1):
        current = extended_results[i]
        next_res = extended_results[i + 1]
        
        # Переход через 0.7
        if load_at_07 is None and current['max_tsai_wu'] < 0.7 <= next_res['max_tsai_wu']:
            # Линейная интерполяция
            t = (0.7 - current['max_tsai_wu']) / (next_res['max_tsai_wu'] - current['max_tsai_wu'])
            load_at_07 = current['load_kg'] + t * (next_res['load_kg'] - current['load_kg'])
        
        # Переход через 0.9
        if load_at_09 is None and current['max_tsai_wu'] < 0.9 <= next_res['max_tsai_wu']:
            t = (0.9 - current['max_tsai_wu']) / (next_res['max_tsai_wu'] - current['max_tsai_wu'])
            load_at_09 = current['load_kg'] + t * (next_res['load_kg'] - current['load_kg'])
        
        # Переход через 1.0
        if load_at_10 is None and current['max_tsai_wu'] < 1.0 <= next_res['max_tsai_wu']:
            t = (1.0 - current['max_tsai_wu']) / (next_res['max_tsai_wu'] - current['max_tsai_wu'])
            load_at_10 = current['load_kg'] + t * (next_res['load_kg'] - current['load_kg'])
    
    print(f"\n  Нагрузка, при которой индекс Цая-Ву достигает 0.7 (зона внимания):")
    print(f"    ≈ {load_at_07:.0f} кг" if load_at_07 else "    не достигнута в рассмотренном диапазоне")
    
    print(f"\n  Нагрузка, при которой индекс Цая-Ву достигает 0.9 (критическое состояние):")
    print(f"    ≈ {load_at_09:.0f} кг" if load_at_09 else "    не достигнута в рассмотренном диапазоне")
    
    print(f"\n  Нагрузка, при которой индекс Цая-Ву достигает 1.0 (прогнозируемое разрушение):")
    print(f"    ≈ {load_at_10:.0f} кг" if load_at_10 else "    не достигнута в рассмотренном диапазоне")
    

    # Результаты прогнозирования
    
    print("\n" + "─" * 80)
    print("ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ")
    print("─" * 80)
    
    loads_all = loads_kg_exp + loads_kg_ext
    sigma_all = [load_to_stress(l) for l in loads_all]
    
    # Для экспериментальных нагрузок используем рассчитанные напряжения
    # Для расширенных — те же самые
    tsai_wu_all = []
    
    # Сначала добавляем результаты экспериментальных нагрузок
    for load in loads_kg_exp:
        sigma = load_to_stress(load)

        F = load * 9.81
        sigma_calc = F * L * y / I
        tsai_wu_all.append(sigma_calc * 0.01) 
    
    for res in extended_results:
        tsai_wu_all.append(res['max_tsai_wu'])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # График 1: Индекс Цая-Ву в зависимости от нагрузки
    ax1.plot(loads_kg_exp, [0.006, 0.024, 0.054, 0.096, 0.121], 'bo-', 
             linewidth=2, markersize=8, label='Экспериментальные данные')
    
    loads_ext_plot = [res['load_kg'] for res in extended_results]
    tsai_wu_ext_plot = [res['max_tsai_wu'] for res in extended_results]
    ax1.plot(loads_ext_plot, tsai_wu_ext_plot, 'ro-', 
             linewidth=2, markersize=6, label='Прогнозируемые значения')
    
    ax1.axhline(y=0.7, color='orange', linestyle='--', linewidth=1.5, label='Зона внимания (0.7)')
    ax1.axhline(y=0.9, color='orange', linestyle=':', linewidth=1.5, label='Критическое состояние (0.9)')
    ax1.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='Разрушение (1.0)')
    
    ax1.set_xlabel('Нагрузка, кг')
    ax1.set_ylabel('Индекс Цая-Ву')
    ax1.set_title('Зависимость индекса Цая-Ву от нагрузки')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # График 2: Напряжения в опасном слое
    sigma_critical = [res['worst_layer_stress'] for res in extended_results]
    ax2.plot(loads_ext_plot, sigma_critical, 'g-s', linewidth=2, markersize=6, 
             label='Напряжение в опасном слое (σ₂)')
    ax2.axhline(y=STRENGTH_TRANSVERSE, color='red', linestyle='--', linewidth=1.5, 
                label=f'Предел прочности поперек ({STRENGTH_TRANSVERSE} МПа)')
    
    ax2.set_xlabel('Нагрузка, кг')
    ax2.set_ylabel('Напряжение, МПа')
    ax2.set_title('Напряжения в наиболее опасном слое (угол 90°)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('extended_prediction.png', dpi=150)
    plt.show()


if __name__ == "__main__":
    main()