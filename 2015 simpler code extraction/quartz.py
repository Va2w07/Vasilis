import extraction
import csv
import numpy as np
import matplotlib.pyplot as plt

def quartz(ref, samp, thickness):
    # create two arrays tArray 0 and 1 with the data of the ref file
    tArray = [list(), list(), list(), list()]
    with open(ref) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            tArray[0].append(float(row[0]))
            tArray[1].append(float(row[1]))

    # create two arrays tArray 2 and 3 with the data of the sample file
    with open(samp) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            tArray[2].append(float(row[0]))
            tArray[3].append(float(row[1]))

    print(tArray[0])
    print(tArray[1])
    print(tArray[2])
    print(tArray[3])

    exp = int(np.log2(len(tArray[0]))) + 5
    L = 2 ** exp
    dt = ((tArray[0][-1] - tArray[0][0]) * 1e-12 / len(tArray[0]))  # average dt
    print(dt)
    xmaxr = np.argmax(tArray[1])
    xmaxs = np.argmax(tArray[3])
    offset = np.abs(int((tArray[0][0] - tArray[2][0]) * 1e-12 / dt))
    tch = ((np.abs(((tArray[0][0] - tArray[2][0]) * 1e-12 / dt))) - (np.abs(int((tArray[0][0] - tArray[2][0]) * 1e-12 / dt))))
    ddt = tch * dt
    print(ddt)
    t1 = np.arange(0, L, 1) * dt
    t2 = np.arange(offset, L + offset, 1) * dt
    wr = window(xmaxr, tArray[1])
    ws = window(xmaxs, tArray[3])
    ref = (np.array(tArray[1]) - tArray[1][0]) * wr
    sample = (np.array(tArray[3]) - tArray[3][0]) * ws

    fref = np.interp(np.arange(0, len(ref)) * dt, (np.array(tArray[0]) - tArray[0][0]) * 1e-12, ref)
    fsam = np.interp(np.arange(0, len(sample)) * dt, (np.array(tArray[2]) - tArray[2][0]) * 1e-12, sample)

    plt.plot(tArray[0], ref, tArray[2], sample)

    ref = fref
    sample = fsam
    (xRVal2, yRVal2, xSVal2, ySVal2, nReal, nImag, f, start, end) = extraction.getRefrac(t1, ref, t2, sample, thickness, 1, ddt, dt, L)
    plt.show()
    return (nReal, nImag, f)

def window(xmax, x):
    width = 150.0
    sigma = width / np.sqrt(2 * np.log(2))
    w = np.ones(len(x))
    for i in range(len(w)):
        if i > (xmax + 20):
            w[i] = np.exp((-(i - xmax) ** 2) / (2 * sigma ** 2))
    return w

n, ni, fx = quartz('/Users/vasilis/Dropbox/Mac/Documents/GitHub/Vasilis/2015 simpler code extraction/2mmr_1.csv', '/Users/vasilis/Dropbox/Mac/Documents/GitHub/Vasilis/2015 simpler code extraction/2mms_1.csv', 2.000)

plt.plot(fx * 1e-12, n, linewidth=2.0, label='real ref index', color='red', linestyle='-')
plt.show()
plt.plot(fx * 1e-12, ni, linewidth=2.0, label='imag ref index', color='blue', linestyle='--')
plt.show()

# nh, nih, fh = quartz('slabhept1.txt', 'hept1.txt', 2.292)
# plt.plot(fh, nh, linewidth=2.0, label='real ref index', color='red', linestyle='-')
# plt.show()
# plt.plot(fh, nih, linewidth=2.0, label='imag ref index', color='blue', linestyle='--')
# plt.show()