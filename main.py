import numpy as np
import requests
import matplotlib.pyplot as plt

from material_dict import *


def get_material_id(name):
    return material_dict[name]


def get_data(name, pressure):
    # pressure unit is MPa
    material_id = get_material_id(name)
    url_base = "http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&Type=IsoBar&Digits=12&THigh=2000&TLow=0&TInc=1&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
    url = url_base + '&ID=' + material_id + '&P=' + str(pressure)

    r = requests.get(url).text.split('\n')

    keys = r[0].split('\t')

    data = []
    for i, e in enumerate(r):
        if i == 0 or e == "":
            continue

        def convert(x):
            if x[0].isdigit() or x[0] == '-':
                return float(x)
            elif x == u'undefined':
                return float('inf')
            else:
                return x
        items = map(convert, e.split('\t'))

        data.append({keys[j]: items[j] for j in range(len(keys))})

    for e in data:
        e[u'log Volume'] = np.log(e[u'Volume (l/mol)'])

    data = calc_thermal_expansion(data)

    return data


def calc_thermal_expansion(data):
    V = u'Volume (l/mol)'
    T = u'Temperature (K)'
    alpha = u'Thermal expansion (/K)'

    for i, e in enumerate(data):
        if i == len(data) - 1:
            e[alpha] = data[i - 1][alpha]
        else:
            try:
                tmp = (data[i + 1][V] - data[i][V]) / (data[i + 1][T] - data[i][T]) / data[i][V]
                e[alpha] = tmp
            except:
                e[alpha] = float('inf')
    return data


def find_T_inverse(data):
    T = u'Temperature (K)'
    JT = u'Joule-Thomson (K/MPa)'
    res = []
    for i, e in enumerate(data):
        if i == len(data) - 1:
            continue
        if data[i][JT] * data[i + 1][JT] < 0.0:
            res.append(data[i][T])
    return res


def find_T_2(data):
    T = u'Temperature (K)'
    alpha = u'Thermal expansion (/K)'
    res = []
    for i, e in enumerate(data):
        if i == len(data) - 1:
            continue
        if (data[i][alpha] * data[i][T] - 1.0) * (data[i + 1][alpha] * data[i][T] - 1.0) < 0.0:
            res.append(data[i][T])
    return res


def get_T1T2(data):
    T1 = find_T_inverse(data)
    T2 = find_T_2(data)
    print T1, T2
    if len(T1) == 2 and len(T2) >= 2:
        if np.abs(T1[0] - T2[0]) < 10:
            return [[T1[0], T2[0]], [T1[1], T2[-1]]]
    elif len(T1) == 1 and len(T2) >= 1:
        if np.abs(T1[0] - T2[0]) < 10:
            return [[T1[0], T2[0]]]
        elif np.abs(T1[0] - T2[-1]) < 10:
            return [[T1[0], T2[-1]]]
    return [[]]


def get_material_data(name):
    d = []
    for p in np.arange(0, 10, 1):
        data = get_data(name, p)
        tmp = get_T1T2(data)
        for e in tmp:
            d.append(e)

    T1 = []
    T2 = []
    for e in d:
        if len(e) < 2:
            continue
        T1.append(e[0])
        T2.append(e[1])
    return T1, T2


def plot_by_pressure(name):
    plt.xlabel('pressure [MPa]')
    plt.ylabel('temperature [K]')
    plt.title(name)

    T1_x = []
    T1_y = []
    T2_x = []
    T2_y = []

    for p in np.arange(0, 10, 1):
        data = get_data(name, p)
        T1 = find_T_inverse(data)
        T2 = find_T_2(data)
        print p, T1, T2
        for e in T1:
            T1_x.append(p)
            T1_y.append(e)
        for e in T2:
            T2_x.append(p)
            T2_y.append(e)

    plt.scatter(T1_x, T1_y, c='red', marker='x', label='T1', s=30)
    plt.scatter(T2_x, T2_y, c='blue', label='T2', s=10)
    plt.legend()
    plt.show()


def parametric_plot():
    T1_N, T2_N = material_dict['Nitrogen']
    T1_W, T2_W = material_dict['Water']
    T1_X, T2_X = material_dict['Xenon']
    T1_T, T2_T = material_dict['Oxygen']

    plt.scatter(T1_N, T2_N, c='red', label='Nitrogen')
    plt.scatter(T1_W, T2_W, c='blue', label='Water')
    plt.scatter(T1_X, T2_X, c='yellow', label='Xenon')
    plt.scatter(T1_T, T2_T, c='green', label='Oxygen')

    plt.xlim([0, 800])
    plt.ylim([0, 800])
    plt.xlabel('T1 [K]')
    plt.ylabel('T2 [K]')
    plt.legend()
    plt.title('1 - 10 MPa')

    x = np.arange(0, 800)
    plt.plot(x, x, 'black')

    plt.show()


if __name__ == '__main__':
    d = get_data('Water', 10)
    get_T1T2(d)
