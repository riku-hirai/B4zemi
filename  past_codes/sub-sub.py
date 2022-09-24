import numpy as np
import csv
import copy

C = 0.5
imax = 30
x = 0.0
dx = 0.1
dt = 0.0
dt_i = 0.0
dt_t = 0.0
tmax = 1
m = 0
mmax = 5000
ary_0 = np.zeros((imax+1,8)) # row array
ary_m = np.zeros((mmax+1,imax+1,8)) # step array


def init(imax):
    for i in range (0, imax+1):

        x = dx * i
        rho = 1.0 - 0.023 * x
        T = 1.0 - 0.009333 * x
        V = 0.05 + 0.11 * x
        A = 0.0
        if i < 16:
            A = 1 + 2.2 * (x - 1.5) ** 2
        else:
            A = 1 + 0.2223 * (x - 1.5) ** 2
        ary_0[i] = [i+1,A,rho,V,T,0,0,0]
    dt_t = calc_dt(ary_0)
    ary_0[0][5] = dt_t


def csv_make(array):
    header = ['グリッドi', '断面積比A', '密度rho', '速度V', '温度T', '圧力P','マッハ数M', '流量m']
    body = array
    
    with open(outputname(mmax), 'w',newline="",encoding='utf_8_sig') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(body)
    f.close()

def outputname(m):
    return 'sub-sub'+ str(m)+ '.csv'

def calc_dt(array):
    for i in range (0, imax+1):
        dt_i = C * dx / (np.sqrt(array[i][4])+array[i][3])
        if i == 0:
            dt_t = dt_i
        else:
            if dt_i < dt_t:
                dt_t = dt_i
    return dt_t

def calculation(m):
    calcA = np.zeros((imax+1,8))
    rowA  = np.zeros((imax+1,8))
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
    dt_step = 0.0
    if m == 1:
        rowA = copy.copy(ary_0)
    else:
        # print(str(ary_m[1]))
        # print(str(ary_m[m-1]))
        for i in range(imax+1):
            for j in range (8):
                rowA[i][j] = copy.copy(ary_m[m-1][i][j])
            
    dt_step = calc_dt(rowA)
    for i in range(0, imax+1):
        calcA[i][0] = copy.copy(rowA[i][0])
        calcA[i][1] = copy.copy(rowA[i][1])
        
        if i < imax:
            rho_pre[i] = -rowA[i][2]*(rowA[i+1][3]-rowA[i][3])/0.1 - rowA[i][2]*rowA[i][3]*(np.log(rowA[i+1][1]) - np.log(rowA[i][1]))/0.1 - rowA[i][3]*(rowA[i+1][2]-rowA[i][2])/0.1 #式合致
            V_pre[i] = -rowA[i][3] * (rowA[i+1][3]-rowA[i][3]) / 0.1 -  1 / 1.4 * ((rowA[i+1][4]-rowA[i][4])/0.1 + rowA[i][4] / rowA[i][2] * (rowA[i+1][2]-rowA[i][2])/0.1)#式合致
            T_pre[i] = - rowA[i][3] * (rowA[i+1][4] - rowA[i][4]) / 0.1 - 0.4 * rowA[i][4] * ((rowA[i+1][3]-rowA[i][3])/0.1 + rowA[i][3]*(np.log(rowA[i+1][1]) - np.log(rowA[i][1]))/0.1) #式合致
            rho_bar[i] = rowA[i][2] + rho_pre[i] * dt_step #式合致
            V_bar[i] = rowA[i][3] + V_pre[i] * dt_step #式合致
            T_bar[i] = rowA[i][4] + T_pre[i] * dt_step #式合致

    for i in range(1, imax):
        if i == 0:
            continue
        else:
            rho_cor[i] = - rho_bar[i] * (V_bar[i]-V_bar[i-1]) / 0.1 - rho_bar[i] * V_bar[i] * (np.log(rowA[i][1]) - np.log(rowA[i-1][1])) / 0.1 -  V_bar[i] * (rho_bar[i] - rho_bar[i-1]) / 0.1
            V_cor[i] = -V_bar[i]*(V_bar[i]-V_bar[i-1])/0.1 - 1/1.4*((T_bar[i]-T_bar[i-1])/0.1+T_bar[i]/rho_bar[i]*(rho_bar[i]-rho_bar[i-1])/0.1) #式合致
            T_cor[i] = -V_bar[i]*(T_bar[i]-T_bar[i-1])/0.1 - 0.4*T_bar[i]*((V_bar[i]-V_bar[i-1])/0.1 + V_bar[i]*(np.log(rowA[i][1]) - np.log(rowA[i-1][1]))/0.1)
            rho_av[i] = 0.5 * (rho_pre[i] + rho_cor[i]) #式合致
            V_av[i] = 0.5 * (V_pre[i] + V_cor[i])
            T_av[i] = 0.5 * (T_pre[i] + T_cor[i])
            calcA[i][2] = rowA[i][2] + rho_av[i] * dt_step
            calcA[i][3] = rowA[i][3] + V_av[i] * dt_step
            calcA[i][4] = rowA[i][4] + T_av[i] * dt_step

    calcA[0][2] = 1
    calcA[0][3] = 2 * calcA[1][3] - calcA[2][3]
    calcA[0][4] = 1

    calcA[imax][2] = 2 * calcA[imax-1][2] - calcA[imax-2][2]
    calcA[imax][3] = 2 * calcA[imax-1][3] - calcA[imax-2][3]
    calcA[imax][4] = 0.93 /calcA[imax][2]

    for i in range(0, imax+1):
        if m ==1:
            if i < imax:
                calcA[i][5] = calcA[i][2] * calcA[i][4]
            else: 
                calcA[i][5] = 0.93
        else:
            if i < imax:
                calcA[i][5] = calcA[i][2] * calcA[i][4]
            else: 
                calcA[i][5] = rowA[i][5]
        calcA[i][6] = calcA[i][3] / np.sqrt(calcA[i][4])
        calcA[i][7] = calcA[i][1] * calcA[i][2] * calcA[i][3]

    for i in range(0, imax+1):
        for j in range(0, 8):
            ary_m[m][i][j] = copy.copy(calcA[i][j])
    # print(str(ary_m[m]))


init(imax)
for m in range(1,mmax+1):
    calculation(m)


csv_make(ary_m[mmax])