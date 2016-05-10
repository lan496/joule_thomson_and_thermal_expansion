import numpy as np
import scipy
import requests
import matplotlib.pyplot as plt

from scipy import optimize
from scipy import interpolate as intp

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


def spline_coefficient(f, x0, x1):
    s = np.linspace(x0, x1, num=4)

    A = np.zeros((4, 4))
    b = np.zeros(4)

    for i in np.arange(4):
        for j in np.arange(4):
            A[i, j] = s[i] ** (3 - j)
        b[i] = f(s[i])

    Ainv = np.linalg.inv(A)

    x = np.dot(Ainv, b)

    return x[0], x[1], x[2], x[3]


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


def divide_data_by_phase(data):
    vaper = filter(lambda e: e['Phase'] == 'vaper', data)
    liquid = filter(lambda e: e['Phase'] == 'liquid', data)
    supercritical = filter(lambda e: e['Phase'] == 'supercritical', data)
    res = []
    for e in [vaper, liquid, supercritical]:
        if len(e) > 0:
            res.append(e)
    return res


def find_log_T1(data):
    res = []

    for i in np.arange(len(data) - 1):
        if data[i]['JT'] * data[i + 1]['JT'] > 0.0:
            continue
        log_T_lst = [e['log_T'] for e in data]
        JT_lst = [e['JT'] for e in data]

        log_T1 = liner_interpolation(log_T_lst[i], log_T_lst[i + 1], JT_lst[i], JT_lst[i + 1])
        err_ip = np.abs(second_order_differential(log_T_lst, JT_lst, i) / 8.0 * ((log_T_lst[i + 1] - log_T_lst[i]) ** 2))

        res.append((log_T1, err_ip))

    return res


def find_log_T2(data):
    x = []
    y = []

    log_T_lst = [e['log_T'] for e in data]
    log_V_lst = [np.log(e['V']) for e in data]

    for i in np.arange(len(data) - 1):
        xi = (log_T_lst[i] + log_T_lst[i + 1]) / 2.0
        try:
            yi = (log_V_lst[i + 1] - log_V_lst[i]) / (log_T_lst[i + 1] - log_T_lst[i]) - 1.0
        except:
            yi = float('inf')
        x.append(xi)
        y.append(yi)

    res = []
    for i in np.arange(len(y) - 1):
        if y[i] * y[i + 1] > 0.0:
            continue

        log_T2 = liner_interpolation(x[i], x[i + 1], y[i], y[i + 1])
        err_diff = np.abs(third_order_differential(log_T_lst, log_V_lst, i) / 24.0 * ((log_T_lst[i + 1] - log_T_lst[i]) ** 2))
        err_ip = np.abs(second_order_differential(x, y, i) / 8.0 * ((x[i + 1] - x[i]) ** 2))

        res.append((log_T2, err_diff + err_ip))

    return res


def find_log_T1_with_spline(data):
    log_T_lst = [e['log_T'] for e in data]
    JT_lst = [e['JT'] for e in data]
    res = []

    f = intp.interp1d(log_T_lst, JT_lst, kind='cubic')

    for i in np.arange(len(log_T_lst) - 1):
        if JT_lst[i] * JT_lst[i + 1] > 0.0:
            continue
        log_T1 = optimize.bisect(f, log_T_lst[i], log_T_lst[i + 1])

        a, b, c, d = spline_coefficient(f, log_T_lst[i], log_T_lst[i + 1])
        err_ip = np.abs(max(6.0 * a * log_T_lst[i] + 2.0 * b, 6.0 * a * log_T_lst[i + 1], 2.0 * b) / 8.0 * ((log_T_lst[i + 1] - log_T_lst[i]) ** 2))

        res.append((log_T1, err_ip))

    return res


def find_log_T2_with_spline(data):
    x = []
    y = []

    log_T_lst = [e['log_T'] for e in data]
    log_V_lst = [np.log(e['V']) for e in data]

    for i in np.arange(len(data) - 1):
        xi = (log_T_lst[i] + log_T_lst[i + 1]) / 2.0
        yi = (log_V_lst[i + 1] - log_V_lst[i]) / (log_T_lst[i + 1] - log_T_lst[i]) - 1.0

        x.append(xi)
        y.append(yi)

    f = intp.interp1d(log_T_lst, log_V_lst, kind='cubic')
    g = intp.interp1d(x, y, kind='cubic')

    res = []

    for i in np.arange(len(x) - 1):
        if y[i] * y[i + 1] > 0.0:
            continue
        log_T2 = optimize.bisect(g, x[i], x[i + 1])
        af, bf, cf, df, = spline_coefficient(f, log_T_lst[i], log_T_lst[i + 1])
        ag, bg, cg, dg = spline_coefficient(g, x[i], x[i + 1])
        err_diff = np.abs(6.0 * af / 24.0 * ((log_T_lst[i + 1] - log_T_lst[i]) ** 2))
        err_ip = np.abs(max(6.0 * ag * x[i] + 2.0 * bg, 6.0 * ag * x[i + 1], 2.0 * bg)) / 8.0 * ((x[i + 1] - x[i]) ** 2)

        res.append((log_T2, err_diff + err_ip))

    return res


def calc_epsilon(log_T1, log_T2):
    if len(log_T1) != len(log_T2):
        return []
    res = []
    for i in np.arange(len(log_T1)):
        T1 = np.exp(log_T1[i][0])
        T2 = np.exp(log_T2[i][0])
        epsilon = np.abs(T1 - T2) / T1
        eps_err = T2 / T1 * max(np.abs(1 - np.exp(log_T1[i][1] + log_T2[i][1])), np.abs(1 - np.exp(-(log_T1[i][1] + log_T2[i][1]))))
        res.append([T1, T2, epsilon, eps_err])
    return res


if __name__ == '__main__':
    print 'Nitrogen'
    for p in np.arange(1, 20, 1):
        data = get_data('Nitrogen', p)
        pdata = divide_data_by_phase(data)
        print str(p) + ' MPa'
        for d in pdata:
            print "naive differential"
            log_T1 = find_log_T1(d)
            log_T2 = find_log_T2(d)
            ep = calc_epsilon(log_T1, log_T2)
            for e in ep:
                print e
            print "third spline"
            log_T1_sp = find_log_T1_with_spline(d)
            log_T2_sp = find_log_T2_with_spline(d)
            ep_sp = calc_epsilon(log_T1_sp, log_T2_sp)
            for e in ep_sp:
                print e
        print
