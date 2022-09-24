import numpy as np
import csv
import copy

x = 0.0
dt = 0.0
dt_t = 0.0
m = 0
gamma = 1.4

C = 0.5
C_x = 0.2
imax = 60 #格子数
jmax = 10 #計算の配列の要素数
dx = 3.0 / imax #長さが3の場合
mmax = 1400 #時間のステップ数
ary_0 = np.zeros((imax+1,jmax)) # row array　初期条件での要素
ary_m = np.zeros((mmax+1,imax+1,jmax)) # step array m回目のステップでの要素
ary_c = np.zeros((mmax+1,imax+1, 7)) # U1, U2, U3, S1, S2, S3, J2
dt_i = 0.0
dt_m = np.zeros(mmax+1)

ary_title = 'shock with no S'
ary_header = ['グリッドi', '断面積比A','U1','U2','U3','密度rho', '速度V', '温度T','マッハ数M', '圧力P']

def init(imax):
    U_1 = 0.0
    U_2 = 0.0
    U_3 = 0.0
    for i in range (0, imax+1):

        x = dx * i
        if i < 10:
            rho = 1.0
            T = 1.0
        elif i < 30:
            rho = 1.0 - 0.366 * (x - 0.5)
            T = 1.0 - 0.167 * (x-0.5)
        elif i < 42:
            rho = 0.634 - 0.702 * (x - 1.5)
            T = 0.833 - 0.4908 * (x - 1.5) 
        else:
            rho = 0.5892 + 0.10228 * (x - 2.1)
            T = 0.93968 + 0.0622 * (x - 2.1)
        A = 1 + 2.2 * (x - 1.5) ** 2
        U_1 = rho * A
        U_2 = 0.59 
        V = U_2 / U_1
        U_3 = U_1 * (T /(gamma  - 1) + gamma / 2 * V**2)
        
        if i < imax:
            P = rho * T
        else:
            P = 0.6784
        ary_0[i] = [i+1,A,U_1,U_2,U_3,rho,V,T,0,P]
    dt_m[m] = calc_dt(ary_0)


def csv_make(strT,arrayH,arrayB):
    #
    title = strT
    header = arrayH
    body = arrayB
    
    with open(outputname(title,mmax), 'w',newline="",encoding='utf_8_sig') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(body)
    f.close()

def outputname(strT,m):
    return strT+ str(m)+ '.csv'

def calc_dt(array):
    for i in range (0, imax+1):
        dt_i = C * dx / (np.sqrt(array[i][7])+array[i][6])
        if i == 0:
            dt_t = dt_i
        else:
            if dt_i < dt_t:
                dt_t = dt_i
    return dt_t

def calc_flux(array):
    for i in range (0, imax+1):
        # F1 F2 F3 S
        ary_c[m][i][0] = array[i][3]
        ary_c[m][i][1] = array[i][3] ** 2 / array[i][2] + (gamma  - 1)/gamma * (array[i][4] - gamma / 2 * array[i][3] ** 2 / array[i][2])
        ary_c[m][i][2] = gamma * array[i][3] * array[i][4] / array[i][2] - gamma * (gamma  - 1) / 2 * array[i][3] ** 3 /array[i][2] ** 2

        if i ==0:
            continue
        elif i < imax:
            ary_c[m][i][3] = C_x * np.abs(array[i+1][8] - 2 * array[i][8] + array[i-1][8]) / (array[i+1][8] + 2 * array[i][8] + array[i-1][8]) * (array[i+1][2] - 2 * array[i][2] + array[i-1][2])
            ary_c[m][i][4] = C_x * np.abs(array[i+1][8] - 2 * array[i][8] + array[i-1][8]) / (array[i+1][8] + 2 * array[i][8] + array[i-1][8]) * (array[i+1][3] - 2 * array[i][3] + array[i-1][3])
            ary_c[m][i][5] = C_x * np.abs(array[i+1][8] - 2 * array[i][8] + array[i-1][8]) / (array[i+1][8] + 2 * array[i][8] + array[i-1][8]) * (array[i+1][4] - 2 * array[i][4] + array[i-1][4])
            ary_c[m][i][6] = 1 / gamma * array[i][5] * array[i][7] * (array[i+1][1] - array[i][1]) / dx

def calc_flux_bar(U1, U2, U3):
    F1 = [0] * (imax+1)
    F2 = [0] * (imax+1)
    F3 = [0] * (imax+1)

    for i in range (0, imax):
            # F1 F2 F3
            F1[i] = U2[i]
            F2[i] = U2[i] ** 2 / U1[i] + (gamma -1)/gamma * (U3[i] - gamma / 2 * U2[i] ** 2 / U1[i])
            F3[i] = gamma * U2[i] * U3[i] / U1[i] - gamma * (gamma  - 1) / 2 * U2[i] ** 3 /U1[i]**2

    return F1,F2,F3

def calculation(m):
    calcA = np.zeros((imax+1,10))
    rowA  = np.zeros((imax+1,10))
    flux_ary = np.zeros((imax+1, 6))
    rho_pre = [0] * (imax+1)
    rho_bar = [0] * (imax+1)
    rho_cor = [0] * (imax+1)
    rho_av = [0] * (imax+1)
    V_pre = [0] * (imax+1)
    V_bar = [0] * (imax+1)
    V_cor = [0] * (imax+1)
    V_av = [0] * (imax+1)
    T_pre = [0] * (imax+1)
    T_bar = [0] * (imax+1)
    T_cor = [0] * (imax+1)
    T_av = [0] * (imax+1)

    rho_shockbar = [0] * (imax+1)
    T_shockbar = [0] * (imax+1)
    F1_bar = [0] * (imax+1)
    F2_bar = [0] * (imax+1)
    F3_bar = [0] * (imax+1)


    if m == 1:
        rowA = copy.copy(ary_0)
    else:
        rowA = copy.copy(ary_m[m-1])

    calc_flux(rowA)        
    dt_m[m] = calc_dt(rowA)

    flux_ary = copy.copy(ary_c[m])



    for i in range(0, imax+1):
        #範囲が0~imaxまで
        #shock でrho V T は　U1 U2 U3
        calcA[i][0] = copy.copy(rowA[i][0])
        calcA[i][1] = copy.copy(rowA[i][1])
        
        if i < imax:
            rho_pre[i] = -(flux_ary[i+1][0] - flux_ary[i][0]) / dx  #
            V_pre[i] = -(flux_ary[i+1][1] - flux_ary[i][1]) / dx  + flux_ary[i][6] #
            T_pre[i] = -(flux_ary[i+1][2] - flux_ary[i][2]) / dx  
            rho_bar[i] = rowA[i][2] + rho_pre[i] * dt_m[m]
            V_bar[i] = rowA[i][3] + V_pre[i] * dt_m[m]
            T_bar[i] = rowA[i][4] + T_pre[i] * dt_m[m] 

            # print(T_bar[i])
            rho_shockbar[i] = rho_bar[i] / rowA[i][1] #
            T_shockbar[i] = (gamma - 1) * (T_bar[i] / rho_bar[i] - gamma / 2 * (V_bar[i] / rho_bar[i]) ** 2) #
    F1_bar, F2_bar, F3_bar = calc_flux_bar(rho_bar, V_bar, T_bar) #

    for i in range(1, imax):
        #範囲が1~imax-1まで
        #shock でrho V T は　U1 U2 U3
        if i == 0:
            continue
        else:
            rho_cor[i] = -(F1_bar[i] - F1_bar[i-1]) / dx #
            V_cor[i] = -(F2_bar[i] - F2_bar[i-1]) / dx  + 1 / gamma * rho_shockbar[i] * T_shockbar[i] * (calcA[i+1][1] - calcA[i][1]) / dx #
            T_cor[i] = -(F3_bar[i] - F3_bar[i-1]) / dx #
            rho_av[i] = 0.5 * (rho_pre[i] + rho_cor[i]) #式合致
            V_av[i] = 0.5 * (V_pre[i] + V_cor[i])
            T_av[i] = 0.5 * (T_pre[i] + T_cor[i])
            calcA[i][2] = rowA[i][2] + rho_av[i] * dt_m[m]
            calcA[i][3] = rowA[i][3] + V_av[i] * dt_m[m]
            calcA[i][4] = rowA[i][4] + T_av[i] * dt_m[m]

    #境界条件
    calcA[0][2] = rowA[0][1] #
    calcA[0][3] = 2 * calcA[1][3] - calcA[2][3] #
    calcA[0][4] = calcA[0][2] * (1.0 / (gamma  - 1) + gamma / 2 * (calcA[0][3] / calcA[0][2])**2) #

    calcA[imax][2] = 2 * calcA[imax-1][2] - calcA[imax-2][2] #
    calcA[imax][3] = 2 * calcA[imax-1][3] - calcA[imax-2][3] #
    calcA[imax][4] = 0.6784 * rowA[imax][1] / (gamma  - 1) + gamma / 2 * calcA[imax][3] ** 2 / calcA[imax][2] #

    for i in range(0, imax+1):
        calcA[i][5] = calcA[i][2] / calcA[i][1]
        calcA[i][6] = calcA[i][3] / calcA[i][2]
        if i < imax:
            calcA[i][7] = (gamma  - 1) * (calcA[i][4]/calcA[i][2] - gamma / 2 * calcA[i][6]**2)
            calcA[i][9] = calcA[i][5] * calcA[i][7]
        else: 
            calcA[i][9] = 0.6784
            calcA[i][7] = 0.6784 / calcA[i][5]
        # print(str(calcA[i][7]))
        calcA[i][8] = calcA[i][6] / np.sqrt(calcA[i][7])

    ary_m[m] = copy.copy(calcA)
    # print(str(ary_m[m]))


init(imax)
for m in range(1,mmax+1):
    calculation(m)


csv_make(ary_title,ary_header,ary_m[mmax])