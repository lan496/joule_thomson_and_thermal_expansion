import numpy as np
import requests
import matplotlib.pyplot as plt

from material_dict import *


def liner_interpolation(x0, x1, y0, y1):
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def second_order_differential(x, y, k):
    d1 = (y[k + 1] - y[k]) / (x[k + 1] - x[k])
    d2 = (y[k] - y[k - 1]) / (x[k] - x[k - 1])
    return (d1 - d2) / (x[k + 1] - x[k])


def third_order_differential(x, y, k):
    d1 = second_order_differential(x, y, k + 1)
    d2 = second_order_differential(x, y, k)
    return (d1 - d2) / (x[k + 1] - x[k])


def get_data(name, pressure):
    # pressure unit is MPa
    material_id = material_dict[name]
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
        tmp_dict = {keys[i]: items[i] for i in range(len(keys))}

        data_i = {}
        data_i[u'T'] = tmp_dict[u'Temperature (K)']
        data_i[u'P'] = tmp_dict[u'Pressure (MPa)']
        data_i[u'V'] = tmp_dict[u'Volume (l/mol)']
        data_i[u'JT'] = tmp_dict[u'Joule-Thomson (K/MPa)']
        data_i[u'Phase'] = tmp_dict[u'Phase']

        data_i[u'log_T'] = np.log(tmp_dict[u'Temperature (K)'])

        data.append(data_i)

    return data


def find_log_T1(data):
    res = []

    for i in np.arange(len(data) - 1):
        if data[i]['JT'] * data[i + 1]['JT'] < 0.0:
            log_V_lst = [e['V'] for e in data]
            JT_lst = [e['JT'] for e in data]

            log_T1 = liner_interpolation(data[i]['log_T'], data[i + 1]['log_T'], JT_lst[i], JT_lst[i + 1])
            err_li = second_order_differential(log_V_lst, JT_lst, i) * (log_V_lst[i + 1] - log_V_lst[i]) * (log_V_lst[i + 1] - log_V_lst[i]) / 8.0

            res.append((log_T1, np.abs(err_li)))

    return res


def find_log_T2(data):
    x = []
    y = []
    for i in np.arange(len(data) - 1):
        xi = (data[i]['log_T'] + data[i + 1]['log_T']) / 2.0
        try:
            yi = (np.log(data[i + 1]['V']) - np.log(data[i]['V'])) / (data[i + 1]['log_T'] - data[i]['log_T']) - 1.0
        except:
            yi = float('inf')
        x.append(xi)
        y.append(yi)

    res = []
    for i in np.arange(len(y) - 1):
        if y[i] * y[i + 1] > 0.0:
            continue

        log_V_lst = [np.log(e['V']) for e in data]

        log_T2 = liner_interpolation(x[i], x[i + 1], y[i], y[i + 1])
        err_diff = third_order_differential(x, log_V_lst, i) * (x[i + 1] - x[i]) * (x[i + 1] - x[i]) / 24.0
        err_li = second_order_differential(x, y, i) * (x[i + 1] - x[i]) * (x[i + 1] - x[i]) / 8.0

        res.append((log_T2, np.abs(err_diff) + np.abs(err_li)))

    return res


def calc_epsilon(log_T1, log_T2):
    if len(log_T1) != len(log_T2):
        return []
    res = []
    for i in np.arange(len(log_T1)):
        T1 = np.exp(log_T1[i][0])
        T2 = np.exp(log_T2[i][0])
        print T1, T2
        epsilon = np.abs(T1 - T2) / T1
        err = T2 / T1 * max(np.abs(1 - np.exp(log_T1[i][1] + log_T2[i][1])), np.abs(1 - np.exp(-(log_T1[i][1] + log_T2[i][1]))))
        res.append((epsilon, err))
    return res


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


if __name__ == '__main__':
    for p in np.arange(1, 30, 1):
        data = get_data('Water', p)
        print str(p) + ' MPa'
        log_T1 = find_log_T1(data)
        log_T2 = find_log_T2(data)
        ep = calc_epsilon(log_T1, log_T2)
        print ep
