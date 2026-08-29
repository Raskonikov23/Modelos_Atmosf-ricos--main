import math
import numpy as np

def DATANH(XH):
    """
    THE DATANH FUNCTION CALCULATES THE HYPERBOLIC TANGENT
    ARC OF AN ARGUMENT X, IN DOUBLE PRECISION.
    """
    if abs(XH) > 1.0:
        raise ValueError("INVALID HYPERBOLIC TANGENT ARC ARGUMENT")
    
    ARGU = (1.0 + XH) / (1.0 - XH)
    DATANH = math.log(ARGU) / 2.0
    return DATANH

def TEMLO(ALTU, C):
    """
    THE FUNCTION TEMLO EVALUATES THE TEMPERATURE AT THE LOCAL REGION.
    """
    ZX = 125.0
    ZO = 90.0
    TO = 188.0
    
    HIGX = ALTU - ZX
    HIGO = ALTU - ZO
    TEMLO = TO

    if HIGO == 0.0:
        return TEMLO
    
    if HIGX > 0.0:
        TEMLO = C[6] + C[3] * math.atan(C[4] * HIGX + C[5] * HIGX * HIGX * HIGX)
    else:
        AUXI = HIGX / HIGO
        TEMLO = C[6] + C[0] * math.atan(C[1] * HIGX + C[2] * HIGX * AUXI * AUXI)
    
    return TEMLO

def IMOWEI(TEXO, ALTU):
    """
    THIS SUBROUTINE USES THE NUMERICALLY INTEGRATED MODEL (JACCHIA, 1977-1) 
    FOR THE MOLECULAR WEIGHT OF THE UPPER ATMOSPHERE.
    """
    # Initialize arrays
    AN = [0.0] * 6
    WT = [0.31111111111, 1.4222222222, 0.53333333333, 1.4222222222, 0.31111111111]
    CI = [28.89122, -2.83071E-02, -6.59924E-03, -3.39574E-04, 6.19256E-05, -1.84796E-06]
    
    AP = [-0.38, 0.0, 0.0, 0.0, 0.0, -0.25]
    WM = [4.0026, 31.9988, 28.0134, 39.948, 15.9994, 1.00797]
    AA = [0.0] * 7
    
    # Constants
    CONV = 2.302585092994
    PI = 3.14159265359
    AVOG = 6.02217E+26
    RGAS = 8.31432
    FRC1 = 5.242E-06
    FRC2 = 0.20955
    FRC3 = 0.78110
    FRC4 = 9.343E-03
    RT = 6356.766
    GO = 9.80665
    RHOI = 3.43E-06
    RMLI = 28.960
    BE = 5.5E-05
    ZO = 90.0
    ZX = 125.0
    TO = 188.0
    
    # Main calculation
    POTE = 0.0045 * (TEXO - TO)
    AUXI = POTE / math.sqrt(1.0 + POTE * POTE)
    TX = TO + 110.5 * DATANH(AUXI)
    AA[0] = 2.0 * (TX - TO) / PI
    AA[3] = 2.0 * (TEXO - TX) / PI
    GX = 1.9 * (TX - TO) / (ZX - ZO)
    AA[1] = GX / AA[0]
    AA[2] = 1.7 * AA[1]
    AA[4] = GX / AA[3]
    AUXI = 28.9 / (TEXO ** 0.25)
    PHIS = (10.0 ** (6.90 + AUXI)) / 2.0E+20
    HND5 = (5.94 + AUXI) * CONV
    AA[5] = BE * AA[4]
    AA[6] = TX
    NC = 5

    STEP = 0.05
    RSLT = 0.0
    ZINI = ZO
    ZEND = 100.0

    if ALTU < ZEND:
        ZEND = ALTU

    ZAUX = (ZEND - ZINI) / 4.0
    AUXI = int(ZAUX / STEP)
    if AUXI <= 0.0:
        AUXI = 1.0
    STEP = ZAUX / AUXI

    while abs(ZEND - ZINI) > 0.0001:
        SUMM = 0.0
        ZAUX = ZINI + 2.0 * STEP
        AUXI = 1.0 + ZAUX / RT
        GRAV = GO / (AUXI * AUXI)
        ZAUX = ZINI

        for I in range(5):
            H = ZAUX - ZO
            WEIG = ((((CI[5] * H + CI[4]) * H + CI[3]) * H + CI[2]) * H + CI[1]) * H + CI[0]
            SUMM = SUMM + GRAV * WEIG / TEMLO(ZAUX, AA) * WT[I]
            ZAUX = ZAUX + STEP

        RSLT = RSLT + STEP * SUMM
        ZINI = ZINI + 4.0 * STEP

    RHOS = RHOI * math.exp(-RSLT / RGAS)
    H = ZEND - ZO
    WEIG = ((((CI[5] * H + CI[4]) * H + CI[3]) * H + CI[2]) * H + CI[1]) * H + CI[0]
    TEMF = TEMLO(ZEND, AA)
    DENU = AVOG * RHOS / RMLI * TO / TEMF
    RHOS = DENU * WEIG
    FAT2 = RHOS / RMLI
    AN[0] = math.log(FRC1 * FAT2)
    AN[1] = math.log(FAT2 * (1.0 + FRC2) - DENU)
    AN[2] = math.log(FRC3 * FAT2)
    AN[3] = math.log(FRC4 * FAT2)
    AN[4] = math.log(2.0 * (DENU - FAT2))
    AN[5] = 0.0

    if ALTU <= 100.0:
        goto_800 = True
    else:
        goto_800 = False

    if not goto_800:
        STEP = 0.05
        ZINI = 100.0
        TEMI = TEMF
        RSLT = 0.0
        ZEND = 140.0
        if ALTU < ZEND:
            ZEND = ALTU

        ZAUX = (ZEND - ZINI) / 4.0
        AUXI = int(ZAUX / STEP)
        if AUXI <= 0.0:
            AUXI = 1.0
        STEP = ZAUX / AUXI

        while abs(ZEND - ZINI) > 0.0001:
            ZAUX = ZINI + 2.0 * STEP
            AUXI = 1.0 + ZAUX / RT
            GRAV = GO / (AUXI * AUXI) / RGAS
            SUMM = 0.0
            ZAUX = ZINI

            for I in range(5):
                SUMM = SUMM + WT[I] * GRAV / TEMLO(ZAUX, AA)
                ZAUX = ZAUX + STEP

            RSLT = RSLT + STEP * SUMM
            ZINI = ZINI + 4.0 * STEP

        TEMF = TEMLO(ZEND, AA)
        AUXI = math.log(TEMI / TEMF)

        for I in range(5):
            AN[I] = AN[I] - RSLT * WM[I] + AUXI * (1 + AP[I])

        if ALTU <= 140.0:
            goto_800 = True

    if not goto_800:
        NC = 6
        ZINI = 140.0
        STEP = 5.0
        RSLT = 0.0
        TEMI = TEMF
        ZEND = 500.0

        while abs(ZEND - ZINI) > 0.0001:
            ZAUX = ZINI + 2.0 * STEP
            AUXI = 1.0 + ZAUX / RT
            GRAV = GO / (AUXI * AUXI) / RGAS
            SUMM = 0.0
            ZAUX = ZINI

            for I in range(5):
                SUMM = SUMM + GRAV * WT[I] / TEMLO(ZAUX, AA)
                ZAUX = ZAUX + STEP

            RSLT = RSLT + STEP * SUMM
            ZINI = ZINI + 4.0 * STEP

        TEMF = TEMLO(ZEND, AA)
        AUXI = math.log(TEMI / TEMF)
        AN[5] = HND5

        for I in range(5):
            AN[I] = AN[I] - RSLT * WM[I] + AUXI * (1 + AP[I])

        if ALTU == 500.0:
            goto_800 = True
        else:
            goto_800 = False

    if not goto_800:
        if ALTU < 500.0:
            TEMI = TEMF
            STEP = -5.0
            ZINI = 500.0
            ZEND = ALTU

            ZAUX = (ZEND - ZINI) / 4.0
            AUXI = int(ZAUX / STEP)
            if AUXI <= 0.0:
                AUXI = 1.0
            STEP = ZAUX / AUXI

            while abs(ZEND - ZINI) > 0.0001:
                AL = [0.0] * 6
                AL[0] = AN[0]
                AL[1] = AN[1] - 0.07 * (1.0 + math.tanh(0.18 * (ZINI - 111.0))) * CONV
                AL[2] = AN[2]
                AL[3] = AN[3]
                AL[4] = AN[4] - 0.24 * math.exp(-0.009 * (ZINI - 97.7) * (ZINI - 97.7)) * CONV

                ZAUX = ZINI + 2.0 * STEP
                AUXI = 1.0 + ZAUX / RT
                GRAV = GO / (AUXI * AUXI) / RGAS
                SUMM = 0.0
                SUMO = 0.0
                ZAUX = ZINI

                for I in range(5):
                    SUMM = SUMM + GRAV * WT[I] / TEMLO(ZAUX, AA)
                    SUMO = SUMO + math.exp(AL[I])
                    ZAUX = ZAUX + STEP

                ZINI = ZINI + 4.0 * STEP
                RSLT_val = STEP * SUMM
                TEMF = TEMLO(ZINI, AA)
                TITF = math.log(TEMI / TEMF)
                TEMI = TEMF

                for I in range(5):
                    AN[I] = AN[I] - RSLT_val * WM[I] + TITF * (1.0 + AP[I])

                SOMA = SUMO / math.exp(AN[5]) * PHIS
                ZAUX = ZINI + 2.0 * STEP
                AUXI = 1.0 + ZAUX / RT
                GRAV = GO / (AUXI * AUXI) / RGAS
                ZAUX = ZINI
                SUMM = 0.0
                SUMO = 0.0
                SUMH = 0.0

                for I in range(5):
                    TEMP = TEMLO(ZAUX, AA)
                    SUMM = SUMM + WT[I] / TEMP
                    AUXI_val = WT[I] / math.sqrt(TEMP)
                    SUMO = SUMO + SOMA * AUXI_val
                    SUMH = SUMH + AUXI_val
                    ZAUX = ZAUX + STEP

                AUXI_val = GRAV * SUMM * STEP * WM[5]
                RSLT_val = AUXI_val - TITF * (1.0 + AP[5])
                AUXI_val = 1000.0 * STEP
                RSLS = SUMO * AUXI_val
                AUXI_val = SUMH * PHIS * AUXI_val
                AN[5] = AN[5] - RSLT_val - RSLS - AUXI_val
        else:
            ZINI = 500.0
            STEP = 2.5
            RSLT = 0.0
            RSLS = 0.0
            TEMI = TEMF
            ZEND = ALTU

            ZAUX = (ZEND - ZINI) / 4.0
            AUXI = int(ZAUX / STEP)
            if AUXI <= 0.0:
                AUXI = 1.0
            STEP = ZAUX / AUXI

            while abs(ZEND - ZINI) > 0.0001:
                ZAUX = ZINI + 2.0 * STEP
                AUXI = 1.0 + ZAUX / RT
                GRAV = GO / (AUXI * AUXI) / RGAS
                SUMM = 0.0
                SUMO = 0.0
                ZAUX = ZINI

                for I in range(5):
                    TEMP = TEMLO(ZAUX, AA)
                    SUMM = SUMM + GRAV / TEMP * WT[I]
                    SUMO = SUMO + WT[I] / math.sqrt(TEMP)
                    ZAUX = ZAUX + STEP

                RSLT = RSLT + STEP * SUMM
                RSLS = RSLS + 1000.0 * STEP * PHIS * SUMO
                ZINI = ZINI + 4.0 * STEP

            TEMF = TEMLO(ZEND, AA)
            AUXI = math.log(TEMI / TEMF)

            for I in range(6):
                AN[I] = AN[I] - RSLT * WM[I] + AUXI * (1.0 + AP[I])

            AN[5] = AN[5] - RSLS

    # Label 800 equivalent
    AN[1] = AN[1] - 0.07 * (1.0 + math.tanh(0.18 * (ZEND - 111.0))) * CONV
    AUXI = ZEND - 97.7
    AN[4] = AN[4] - 0.24 * math.exp(-0.009 * AUXI * AUXI) * CONV
    WEIG = 0.0
    SUMM = 0.0

    for I in range(NC):
        AUXI = math.exp(AN[I])
        AN[I] = AN[I] / CONV
        if AN[I] < 0.0:
            AN[I] = 0.0
        WEIG = WEIG + AUXI * WM[I]
        SUMM = SUMM + AUXI

    WMOL = WEIG / SUMM
    RHOD = WEIG / AVOG

    return AN, WMOL, RHOD, TEMF

# Exemplo de uso:
if __name__ == "__main__":
    TEXO = 1000.0  # Temperatura exosférica em Kelvin
    ALTU = 200.0   # Altitude em km
    
    AN, WMOL, RHOD, TEMF = IMOWEI(TEXO, ALTU)
    
    print("Resultados:")
    print(f"AN: {AN}")
    print(f"WMOL: {WMOL}")
    print(f"RHOD: {RHOD}")
    print(f"TEMF: {TEMF}")