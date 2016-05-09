import numpy as np
import requests
import matplotlib.pyplot as plt

from material_dict import *


def liner_interpolation(x0, x1, y0, y1):
    return x0 - y0 * (x1 - x0) / (y1 - y0)


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


def find_T1(data):
    res = []
    for i in np.arange(len(data) - 1):
        if data[i]['JT'] * data[i + 1]['JT'] < 0.0:
            log_T1 = liner_interpolation(data[i]['log_T'], data[i + 1]['log_T'], data[i]['JT'], data[i + 1]['JT'])
            err = 0.0
            res.append((np.exp(log_T1), err))
    return res


def find_T2(data):
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

    plt.plot(x, y)
    plt.pause(0.01)
    res = []
    for i in np.arange(len(y) - 1):
        if y[i] * y[i + 1] > 0.0:
            continue
        log_T2 = liner_interpolation(x[i], x[i + 1], y[i], y[i + 1])
        err = 0.0
        res.append((np.exp(log_T2), err))

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
    plt.xlabel('log V')
    plt.ylabel(r'$\frac{d}{dx}$ log V - 1')
    plt.xlim(0, 8)
    plt.ylim(-0.1, 0.1)
    plt.title('T2')
    x = np.arange(0, 10)
    y = [0 for e in x]
    plt.plot(x, y, color='black')

    for p in np.arange(1, 10, 3):
        data = get_data('Nitrogen', p)
        print p
        print find_T1(data)
        print find_T2(data)
    plt.show()
