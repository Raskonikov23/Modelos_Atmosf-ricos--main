import math
import numpy as np
from IMOWEI import IMOWEI

def calcular_arrasto(x, area, massa, cd=2.2):
    """
    Calcula a aceleração devido ao arrasto atmosférico
    x: vetor de estado [x, y, z, vx, vy, vz]
    area: área de referência (m²)
    massa: massa do satélite (kg)
    cd: coeficiente de arrasto
    """
    # Constantes
    re = 6378165.0  # Raio da Terra em metros
    omega_terra = 7.292115e-5  # Velocidade angular da Terra (rad/s)
    
    # Calcular altitude
    pos = math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    altitude = (pos - re) / 1000.0  # km
    
    # Obter parâmetros atmosféricos do IMOWEI
    # Temperatura exosférica fixa para simplificação (pode ser variável)
    texo = 1000.0  # K
    
    try:
        AN, WMOL, RHOD, TEMF = IMOWEI(texo, altitude)
        
        # Densidade atmosférica em kg/m³
        densidade = RHOD
        
        # Velocidade relativa à atmosfera (considerando rotação da Terra)
        v_rel = [
            x[3] + omega_terra * x[1],
            x[4] - omega_terra * x[0],
            x[5]
        ]
        
        v_mag = math.sqrt(v_rel[0]**2 + v_rel[1]**2 + v_rel[2]**2)
        
        if v_mag < 1e-10:
            return [0.0, 0.0, 0.0]
        
        # Aceleração devido ao arrasto
        fator = -0.5 * cd * area * densidade * v_mag / massa
        
        a_arrasto = [
            fator * v_rel[0],
            fator * v_rel[1],
            fator * v_rel[2]
        ]
        
        return a_arrasto
        
    except Exception as e:
        print(f"Erro no cálculo do arrasto em {altitude:.1f} km: {e}")
        return [0.0, 0.0, 0.0]

def f_orbital_com_arrasto(t, n, x, area=1.0, massa=1000.0, cd=2.2):
    """
    Equações de movimento orbital COM arrasto atmosférico
    """
    mi = 3.9860320e14
    pos = math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    
    xp = [0.0] * n
    xp[0] = x[3]  # dx/dt = vx
    xp[1] = x[4]  # dy/dt = vy
    xp[2] = x[5]  # dz/dt = vz
    
    # Acelerações gravitacionais
    xp[3] = -mi * x[0] / (pos**3)
    xp[4] = -mi * x[1] / (pos**3)
    xp[5] = -mi * x[2] / (pos**3)
    
    # Adicionar arrasto atmosférico (apenas em baixas altitudes)
    altitude = (pos - 6378165.0) / 1000.0
    if altitude < 1000.0:  # Aplicar arrasto apenas abaixo de 1000 km
        a_arrasto = calcular_arrasto(x, area, massa, cd)
        xp[3] += a_arrasto[0]
        xp[4] += a_arrasto[1]
        xp[5] += a_arrasto[2]
    
    return xp