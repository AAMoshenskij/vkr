import math
import matplotlib.pyplot as plt

class Matrix:
    def __init__(self):
        self.data = [[0.0 for _ in range(3)] for _ in range(3)]
    
    def __str__(self):
        """Возвращает строковое представление матрицы"""
        result = ""
        for i in range(3):
            result += "| "
            for j in range(3):
                result += f"{self.data[i][j]:10.4f} "
            result += "|\n"
        return result

def add_matrices(matrices):
    result = Matrix()
    for m in matrices:
        for i in range(3):
            for j in range(3):
                result.data[i][j] += m.data[i][j]
    return result

def multiply_scalar(m, scalar):
    result = Matrix()
    for i in range(3):
        for j in range(3):
            result.data[i][j] = m.data[i][j] * scalar
    return result

def inverse(m):
    result = Matrix()
    det = (m.data[0][0] * (m.data[1][1] * m.data[2][2] - m.data[1][2] * m.data[2][1]) -
           m.data[0][1] * (m.data[1][0] * m.data[2][2] - m.data[1][2] * m.data[2][0]) +
           m.data[0][2] * (m.data[1][0] * m.data[2][1] - m.data[1][1] * m.data[2][0]))

    if abs(det) < 1e-10:
        print("Error: Matrix is not invertible.")
        return result

    inv_det = 1.0 / det

    result.data[0][0] = (m.data[1][1] * m.data[2][2] - m.data[1][2] * m.data[2][1]) * inv_det
    result.data[0][1] = -(m.data[0][1] * m.data[2][2] - m.data[0][2] * m.data[2][1]) * inv_det
    result.data[0][2] = (m.data[0][1] * m.data[1][2] - m.data[0][2] * m.data[1][1]) * inv_det

    result.data[1][0] = -(m.data[1][0] * m.data[2][2] - m.data[1][2] * m.data[2][0]) * inv_det
    result.data[1][1] = (m.data[0][0] * m.data[2][2] - m.data[0][2] * m.data[2][0]) * inv_det
    result.data[1][2] = -(m.data[0][0] * m.data[1][2] - m.data[0][2] * m.data[1][0]) * inv_det

    result.data[2][0] = (m.data[1][0] * m.data[2][1] - m.data[1][1] * m.data[2][0]) * inv_det
    result.data[2][1] = -(m.data[0][0] * m.data[2][1] - m.data[0][1] * m.data[2][0]) * inv_det
    result.data[2][2] = (m.data[0][0] * m.data[1][1] - m.data[0][1] * m.data[1][0]) * inv_det

    return result

def determinant(m):
    det = (m.data[0][0] * m.data[1][1] * m.data[2][2] + 
           m.data[2][0] * m.data[0][1] * m.data[1][2] + 
           m.data[0][2] * m.data[1][0] * m.data[2][1] - 
           (m.data[0][2] * m.data[1][1] * m.data[2][0] + 
            m.data[0][0] * m.data[2][1] * m.data[1][2] + 
            m.data[2][2] * m.data[1][0] * m.data[0][1]))
    return det

def print_matrix(m):
    for i in range(3):
        for j in range(3):
            print(f"{m.data[i][j]:10.4f}", end=" ")
        print()
    print()

def sum_vector(vec):
    total = 0.0
    for value in vec:
        total += value
    return total

def print_vector(vec):
    for element in vec:
        print(element, end=" ")
    print()

def main():
    n = int(input("Введите количество слоев композита: "))
    h = 1.0 / n

    if n == 1:
        plot_choice = input("Будем ли строить график? (y/n): ")
        if plot_choice == 'y':
            E1 = float(input("Введите E1: "))
            E2 = float(input("Введите E2: "))
            G12 = float(input("Введите G12: "))
            u12 = float(input("Введите коэффициент Пуассона u12: "))
            u21 = float(input("Введите коэффициент Пуассона u21: "))

            angles = []
            Ea_values = []
            Eb_values = []
            Gab_values = []
            u_ab_values = []
            
            for angle_deg in range(0, 91):
                angle_rad = angle_deg * math.pi / 180.0
                angles.append(angle_deg)

                G0 = Matrix()
                G0.data[0][0] = E1 / (1 - u12 * u21)
                G0.data[0][1] = u21 * E1 / (1 - u12 * u21)
                G0.data[1][0] = G0.data[0][1]
                G0.data[1][1] = E2 / (1 - u12 * u21)
                G0.data[2][2] = G12
                G0.data[0][2] = G0.data[1][2] = G0.data[2][0] = G0.data[2][1] = 0

                Gxy = Matrix()
                Gxy.data[0][0] = (math.cos(angle_rad)**4 * G0.data[0][0] + 
                                  math.sin(angle_rad)**4 * G0.data[1][1] + 
                                  2 * math.cos(angle_rad)**2 * math.sin(angle_rad)**2 * (G0.data[0][1] + 2 * G0.data[2][2]))
                Gxy.data[0][1] = (math.sin(angle_rad)**2 * math.cos(angle_rad)**2 * (G0.data[0][0] + G0.data[1][1] - 4 * G0.data[2][2]) + 
                                  G0.data[0][1] * (math.cos(angle_rad)**4 + math.sin(angle_rad)**4))
                Gxy.data[0][2] = (math.sin(angle_rad) * math.cos(angle_rad) * 
                                  (math.cos(angle_rad)**2 * G0.data[0][0] - 
                                   math.sin(angle_rad)**2 * G0.data[1][1] + 
                                   (G0.data[0][1] + 2 * G0.data[2][2]) * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)))
                Gxy.data[1][0] = Gxy.data[0][1]
                Gxy.data[1][1] = (math.sin(angle_rad)**4 * G0.data[0][0] + 
                                  math.cos(angle_rad)**4 * G0.data[1][1] + 
                                  2 * math.cos(angle_rad)**2 * math.sin(angle_rad)**2 * (G0.data[0][1] + 2 * G0.data[2][2]))
                Gxy.data[1][2] = (math.sin(angle_rad) * math.cos(angle_rad) * 
                                  (math.sin(angle_rad)**2 * G0.data[0][0] - 
                                   math.cos(angle_rad)**2 * G0.data[1][1] - 
                                   (G0.data[0][1] + 2 * G0.data[2][2]) * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)))
                Gxy.data[2][0] = Gxy.data[0][2]
                Gxy.data[2][1] = Gxy.data[1][2]
                Gxy.data[2][2] = (math.sin(angle_rad)**2 * math.cos(angle_rad)**2 * (G0.data[0][0] - 2 * G0.data[0][1] + G0.data[1][1]) + 
                                  G0.data[2][2] * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)**2)
                
                detGxy = determinant(Gxy)
                Ea = detGxy / (Gxy.data[1][1] * Gxy.data[2][2] - Gxy.data[1][2]**2)
                Eb = detGxy / (Gxy.data[0][0] * Gxy.data[2][2] - Gxy.data[0][2]**2)
                Gab = detGxy / (Gxy.data[0][0] * Gxy.data[1][1] - Gxy.data[0][1]**2)
                u_ab = (Gxy.data[0][1] * Gxy.data[2][2] - Gxy.data[0][2] * Gxy.data[1][2]) / (Gxy.data[1][1] * Gxy.data[2][2] - Gxy.data[1][2]**2)
                Ea_values.append(Ea)
                Eb_values.append(Eb)
                Gab_values.append(Gab)
                u_ab_values.append(u_ab)

            plt.plot(angles, Ea_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль упругости $E_a$ (МПа)")
            plt.show()

            plt.plot(angles, Eb_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль упругости $E_b$ (МПа)")
            plt.show()

            plt.plot(angles, Gab_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль сдвига $G_{ab}$ (МПа)")
            plt.show()

            plt.plot(angles, u_ab_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Коэффициент Пуассона $u_{ab}$")
            plt.show()
            return
    
    if n == 4:
        plot_choice = input("Будем ли строить график? (y/n): ")
        if plot_choice == 'y':
            E1_values = []
            E2_values = []
            for i in range(n):
                E1 = float(input(f"Введите E1 для слоя {i + 1}: "))
                E2 = float(input(f"Введите E2 для слоя {i + 1}: "))
                E1_values.append(E1)
                E2_values.append(E2)

            G12_values = []
            for i in range(n):
                G12 = float(input(f"Введите G12 для слоя {i + 1}: "))
                G12_values.append(G12)

            u12 = float(input("Введите коэффициент Пуассона u12: "))
            u21 = float(input("Введите коэффициент Пуассона u21: "))

            G0_values = []
            for i in range(n):
                G0 = Matrix()
                G0.data[0][0] = E1_values[i] / (1 - u12 * u21)
                G0.data[0][1] = u21 * E1_values[i] / (1 - u12 * u21)
                G0.data[1][0] = G0.data[0][1]
                G0.data[1][1] = E2_values[i] / (1 - u12 * u21)
                G0.data[2][2] = G12_values[i]
                G0.data[0][2] = G0.data[1][2] = G0.data[2][0] = G0.data[2][1] = 0
                G0_values.append(G0)

            angles = []
            Ea_values = []
            Eb_values = []
            Gab_values = []
            u_ab_values = []
            
            for angle_deg in range(0, 91):
                angle_rad = angle_deg * math.pi / 180.0
                angles.append(angle_deg)

                Gxy_values = []
                for i in range(n):
                    Gxy = Matrix()
                    Gxy.data[0][0] = (math.cos(angle_rad)**4 * G0_values[i].data[0][0] + 
                                      math.sin(angle_rad)**4 * G0_values[i].data[1][1] + 
                                      2 * math.cos(angle_rad)**2 * math.sin(angle_rad)**2 * (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]))
                    Gxy.data[0][1] = (math.sin(angle_rad)**2 * math.cos(angle_rad)**2 * (G0_values[i].data[0][0] + G0_values[i].data[1][1] - 4 * G0_values[i].data[2][2]) + 
                                      G0_values[i].data[0][1] * (math.cos(angle_rad)**4 + math.sin(angle_rad)**4))
                    Gxy.data[0][2] = (math.sin(angle_rad) * math.cos(angle_rad) * 
                                      (math.cos(angle_rad)**2 * G0_values[i].data[0][0] - 
                                       math.sin(angle_rad)**2 * G0_values[i].data[1][1] + 
                                       (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]) * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)))
                    Gxy.data[1][0] = Gxy.data[0][1]
                    Gxy.data[1][1] = (math.sin(angle_rad)**4 * G0_values[i].data[0][0] + 
                                      math.cos(angle_rad)**4 * G0_values[i].data[1][1] + 
                                      2 * math.cos(angle_rad)**2 * math.sin(angle_rad)**2 * (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]))
                    Gxy.data[1][2] = (math.sin(angle_rad) * math.cos(angle_rad) * 
                                      (math.sin(angle_rad)**2 * G0_values[i].data[0][0] - 
                                       math.cos(angle_rad)**2 * G0_values[i].data[1][1] - 
                                       (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]) * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)))
                    Gxy.data[2][0] = Gxy.data[0][2]
                    Gxy.data[2][1] = Gxy.data[1][2]
                    Gxy.data[2][2] = (math.sin(angle_rad)**2 * math.cos(angle_rad)**2 * (G0_values[i].data[0][0] - 2 * G0_values[i].data[0][1] + G0_values[i].data[1][1]) + 
                                      G0_values[i].data[2][2] * (math.sin(angle_rad)**2 - math.cos(angle_rad)**2)**2)
                    Gxy_values.append(Gxy)
                
                Gxy_summ = add_matrices(Gxy_values)
                Gxy_summ = multiply_scalar(Gxy_summ, h)
                detGxy = determinant(Gxy_summ)
                Ea = detGxy / (Gxy_summ.data[1][1] * Gxy_summ.data[2][2] - Gxy_summ.data[1][2]**2)
                Eb = detGxy / (Gxy_summ.data[0][0] * Gxy_summ.data[2][2] - Gxy_summ.data[0][2]**2)
                Gab = detGxy / (Gxy_summ.data[0][0] * Gxy_summ.data[1][1] - Gxy_summ.data[0][1]**2)
                u_ab = (Gxy_summ.data[0][1] * Gxy_summ.data[2][2] - Gxy_summ.data[0][2] * Gxy_summ.data[1][2]) / (Gxy_summ.data[1][1] * Gxy_summ.data[2][2] - Gxy_summ.data[1][2]**2)
                Ea_values.append(Ea)
                Eb_values.append(Eb)
                Gab_values.append(Gab)
                u_ab_values.append(u_ab)

            plt.plot(angles, Ea_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль упругости $E_a$ (МПа)")
            plt.show()

            plt.plot(angles, Eb_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль упругости $E_b$ (МПа)")
            plt.show()

            plt.plot(angles, Gab_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Модуль сдвига $G_{ab}$ (МПа)")
            plt.show()

            plt.plot(angles, u_ab_values)
            plt.xlabel("Угол армирования (°)")
            plt.ylabel("Коэффициент Пуассона $u_{ab}$")
            plt.show()
            return

    E1_values = []
    E2_values = []
    for i in range(n):
        E1 = float(input(f"Введите E1 для слоя {i + 1}: "))
        E2 = float(input(f"Введите E2 для слоя {i + 1}: "))
        E1_values.append(E1)
        E2_values.append(E2)

    G12_values = []
    for i in range(n):
        G12 = float(input(f"Введите G12 для слоя {i + 1}: "))
        G12_values.append(G12)

    u12 = float(input("Введите коэффициент Пуассона u12: "))
    u21 = float(input("Введите коэффициент Пуассона u21: "))

    angle_values = []
    for i in range(n):
        angle_deg = float(input(f"Введите угол армирования (в градусах) для слоя {i + 1}: "))
        angle_rad = angle_deg * math.pi / 180.0
        angle_values.append(angle_rad)

    sigma1_values = []
    sigma2_values = []
    tau12_values = []
    for i in range(n):
        sigma1 = float(input(f"Введите sigma1 для слоя {i + 1}: "))
        sigma2 = float(input(f"Введите sigma2 для слоя {i + 1}: "))
        tau12 = float(input(f"Введите tau12 для слоя {i + 1}: "))
        sigma1_values.append(sigma1)
        sigma2_values.append(sigma2)
        tau12_values.append(tau12)

    G0_values = []
    for i in range(n):
        G0 = Matrix()
        G0.data[0][0] = E1_values[i] / (1 - u12 * u21)
        G0.data[0][1] = u21 * E1_values[i] / (1 - u12 * u21)
        G0.data[1][0] = G0.data[0][1]
        G0.data[1][1] = E2_values[i] / (1 - u12 * u21)
        G0.data[2][2] = G12_values[i]
        G0.data[0][2] = G0.data[1][2] = G0.data[2][0] = G0.data[2][1] = 0
        G0_values.append(G0)
        #print(G0)

    Gxy_values = []
    for i in range(n):
        Gxy = Matrix()
        Gxy.data[0][0] = (math.cos(angle_values[i])**4 * G0_values[i].data[0][0] + 
                          math.sin(angle_values[i])**4 * G0_values[i].data[1][1] + 
                          2 * math.cos(angle_values[i])**2 * math.sin(angle_values[i])**2 * (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]))
        Gxy.data[0][1] = (math.sin(angle_values[i])**2 * math.cos(angle_values[i])**2 * (G0_values[i].data[0][0] + G0_values[i].data[1][1] - 4 * G0_values[i].data[2][2]) + 
                          G0_values[i].data[0][1] * (math.cos(angle_values[i])**4 + math.sin(angle_values[i])**4))
        Gxy.data[0][2] = (math.sin(angle_values[i]) * math.cos(angle_values[i]) * 
                          (math.cos(angle_values[i])**2 * G0_values[i].data[0][0] - 
                           math.sin(angle_values[i])**2 * G0_values[i].data[1][1] + 
                           (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]) * (math.sin(angle_values[i])**2 - math.cos(angle_values[i])**2)))
        Gxy.data[1][0] = Gxy.data[0][1]
        Gxy.data[1][1] = (math.sin(angle_values[i])**4 * G0_values[i].data[0][0] + 
                          math.cos(angle_values[i])**4 * G0_values[i].data[1][1] + 
                          2 * math.cos(angle_values[i])**2 * math.sin(angle_values[i])**2 * (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]))
        Gxy.data[1][2] = (math.sin(angle_values[i]) * math.cos(angle_values[i]) * 
                          (math.sin(angle_values[i])**2 * G0_values[i].data[0][0] - 
                           math.cos(angle_values[i])**2 * G0_values[i].data[1][1] - 
                           (G0_values[i].data[0][1] + 2 * G0_values[i].data[2][2]) * (math.sin(angle_values[i])**2 - math.cos(angle_values[i])**2)))
        Gxy.data[2][0] = Gxy.data[0][2]
        Gxy.data[2][1] = Gxy.data[1][2]
        Gxy.data[2][2] = (math.sin(angle_values[i])**2 * math.cos(angle_values[i])**2 * (G0_values[i].data[0][0] - 2 * G0_values[i].data[0][1] + G0_values[i].data[1][1]) + 
                          G0_values[i].data[2][2] * (math.sin(angle_values[i])**2 - math.cos(angle_values[i])**2)**2)
        Gxy_values.append(Gxy)
        #print(Gxy)
        #tan4 = 4 * (Gxy.data[0][2] - Gxy.data[1][2]) / (Gxy.data[0][0] + Gxy.data[1][1] - 2 * Gxy.data[0][1] - 4 * Gxy.data[2][2])
        #tan2 = 2 * (Gxy.data[0][2] + Gxy.data[1][2]) / (Gxy.data[0][0] - Gxy.data[1][1])
        #print(tan4)
        #print(tan2)

    sigmaX_values = []
    sigmaY_values = []
    tauXY_values = []
    for i in range(n):
        sigmaX = (math.cos(angle_values[i])**2 * sigma1_values[i] + 
                  math.sin(angle_values[i])**2 * sigma2_values[i] - 
                  2 * math.sin(angle_values[i]) * math.cos(angle_values[i]) * tau12_values[i])
        sigmaY = (math.sin(angle_values[i])**2 * sigma1_values[i] + 
                  math.cos(angle_values[i])**2 * sigma2_values[i] + 
                  2 * math.sin(angle_values[i]) * math.cos(angle_values[i]) * tau12_values[i])
        tauXY = (math.sin(angle_values[i]) * math.cos(angle_values[i]) * sigma1_values[i] - 
                 math.sin(angle_values[i]) * math.cos(angle_values[i]) * sigma2_values[i] + 
                 (math.cos(angle_values[i])**2 - math.sin(angle_values[i])**2) * tau12_values[i])
        sigmaX_values.append(sigmaX)
        sigmaY_values.append(sigmaY)
        tauXY_values.append(tauXY)

    sigmaX = sum_vector(sigmaX_values) * h
    sigmaY = sum_vector(sigmaY_values) * h
    tauXY = sum_vector(tauXY_values) * h

    Gxy_summ = add_matrices(Gxy_values)
    Sxy_summ = inverse(Gxy_summ)

    Gxy_summ = multiply_scalar(Gxy_summ, h)
    Sxy_summ = multiply_scalar(Sxy_summ, h)
    detGxy = determinant(Gxy_summ)

    Ea = detGxy / (Gxy_summ.data[1][1] * Gxy_summ.data[2][2] - Gxy_summ.data[1][2]**2)
    Eb = detGxy / (Gxy_summ.data[0][0] * Gxy_summ.data[2][2] - Gxy_summ.data[0][2]**2)
    Gab = detGxy / (Gxy_summ.data[0][0] * Gxy_summ.data[1][1] - Gxy_summ.data[0][1]**2)
    u_ab = (Gxy_summ.data[0][1] * Gxy_summ.data[2][2] - Gxy_summ.data[0][2] * Gxy_summ.data[1][2]) / (Gxy_summ.data[1][1] * Gxy_summ.data[2][2] - Gxy_summ.data[1][2]**2)

    eps_x = (Sxy_summ.data[0][0] * sigmaX + Sxy_summ.data[0][1] * sigmaY + Sxy_summ.data[0][2] * tauXY)
    eps_y = (Sxy_summ.data[1][0] * sigmaX + Sxy_summ.data[1][1] * sigmaY + Sxy_summ.data[1][2] * tauXY)
    gamma_xy = (Sxy_summ.data[2][0] * sigmaX + Sxy_summ.data[2][1] * sigmaY + Sxy_summ.data[2][2] * tauXY)

    print(f"Модуль Юнга Ea итогового материала: {Ea:.4f}")
    print(f"Модуль Юнга Eb итогового материала: {Eb:.4f}")
    print(f"Модуль сдвига Gab итогового материала: {Gab:.4f}")
    print(f"Коэффициент Пуассона итогового материала: {u_ab:.4f}")
    print(f"Напряжение sigmaX итогового материала: {sigmaX:.4f}")
    print(f"Напряжение sigmaY итогового материала: {sigmaY:.4f}")
    print(f"Напряжение tauXY итогового материала: {tauXY:.4f}")
    print(f"Деформация eps_x итогового материала: {eps_x:.4f}")
    print(f"Деформация eps_y итогового материала: {eps_y:.4f}")
    print(f"Деформация gamma_xy итогового материала: {gamma_xy:.4f}")

if __name__ == "__main__":
    main()