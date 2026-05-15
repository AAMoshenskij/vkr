"""
Обертка для вызова stress_strain.py
"""
import sys
import io
from contextlib import redirect_stdout, redirect_stderr
import re

def calculate_laminate_properties(E1_list, E2_list, G12_list, nu12, nu21, 
                                   angles_deg, sigma1_list, sigma2_list, tau12_list,
                                   debug=False):
    """
    Вызов исходной программы с заданными параметрами
    """
    
    n = len(E1_list)
    
    # Формируем входные данные для программы
    input_data = []
    
    # Количество слоев
    input_data.append(str(n))
    
    # Для каждого слоя: E1, E2
    for i in range(n):
        input_data.append(str(E1_list[i]))
        input_data.append(str(E2_list[i]))
    
    # Для каждого слоя: G12
    for i in range(n):
        input_data.append(str(G12_list[i]))
    
    # Коэффициенты Пуассона
    input_data.append(str(nu12))
    input_data.append(str(nu21))
    
    # Углы армирования
    for i in range(n):
        input_data.append(str(angles_deg[i]))
    
    # Напряжения в слоях: sigma1, sigma2, tau12
    for i in range(n):
        input_data.append(str(sigma1_list[i]))
        input_data.append(str(sigma2_list[i]))
        input_data.append(str(tau12_list[i]))
    
    input_text = "\n".join(input_data) + "\n"
    
    captured_output = io.StringIO()
    
    original_stdin = sys.stdin
    original_stdout = sys.stdout
    
    try:
        sys.stdin = io.StringIO(input_text)
        sys.stdout = captured_output
        
        import stress_strain
        
        stress_strain.main()
        
    except Exception as e:
    
        return {
            'Ea': None, 'Eb': None, 'Gab': None, 'nu_ab': None,
            'sigmaX': None, 'sigmaY': None, 'tauXY': None,
            'eps_x': None, 'eps_y': None, 'gamma_xy': None
        }
        
    finally:
        # Восстанавливаем stdin/stdout
        sys.stdin = original_stdin
        sys.stdout = original_stdout
    
    output_text = captured_output.getvalue()
    
    results = _parse_output(output_text, debug)

    
    return results


def _parse_output(output_text, debug=False):
    """
    Парсинг вывода программы для извлечения результатов
    """
    results = {
        'Ea': None, 'Eb': None, 'Gab': None, 'nu_ab': None,
        'sigmaX': None, 'sigmaY': None, 'tauXY': None,
        'eps_x': None, 'eps_y': None, 'gamma_xy': None
    }
    
    # Ищем числовые значения в выводе
    patterns = {
        'Ea': r"Модуль Юнга Ea итогового материала:\s*([-\d.]+)",
        'Eb': r"Модуль Юнга Eb итогового материала:\s*([-\d.]+)",
        'Gab': r"Модуль сдвига Gab итогового материала:\s*([-\d.]+)",
        'nu_ab': r"Коэффициент Пуассона итогового материала:\s*([-\d.]+)",
        'sigmaX': r"Напряжение sigmaX итогового материала:\s*([-\d.]+)",
        'sigmaY': r"Напряжение sigmaY итогового материала:\s*([-\d.]+)",
        'tauXY': r"Напряжение tauXY итогового материала:\s*([-\d.]+)",
        'eps_x': r"Деформация eps_x итогового материала:\s*([-\d.]+)",
        'eps_y': r"Деформация eps_y итогового материала:\s*([-\d.]+)",
        'gamma_xy': r"Деформация gamma_xy итогового материала:\s*([-\d.]+)"
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, output_text)
        if match:
            try:
                results[key] = float(match.group(1))
                if debug:
                    print(f"  Найдено {key}: {match.group(1)}")
            except ValueError:
                if debug:
                    print(f"  Ошибка преобразования {key}: {match.group(1)}")
        else:
            if debug:
                print(f"  Не найдено {key}")
    
    return results

if __name__ == "__main__":
    # Тестовые данные
    n = 2
    E1_list = [46000, 46000]
    E2_list = [18000, 18000]
    G12_list = [4500, 4500]
    nu12 = 0.2
    nu21 = 0.2
    angles_deg = [90, 90]
    sigma1_list = [100, 100]
    sigma2_list = [50, 50]
    tau12_list = [20, 20]
    
    results = calculate_laminate_properties(
        E1_list, E2_list, G12_list, nu12, nu21,
        angles_deg, sigma1_list, sigma2_list, tau12_list,
        debug=True
    )
    
    print("\nРезультаты из обертки:")
    for key, value in results.items():
        print(f"  {key}: {value}")