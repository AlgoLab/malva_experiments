import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
rc('text')
font = {'family' : 'normal',
        'size'   : 17}
matplotlib.rc('font', **font)

def parse_csv(fpath):
    times = []
    ram = []
    y_labels = []
    for line in open(fpath, 'r'):
        if line[0] == 'X':
            continue
        line = line.strip("\n").split(',')
        if line[0] == "RAM":
            ram = [int(x) for x in line[1:]]
        else:
            y_labels.append(line[0])
            times.append([int(x)/60 for x in line[1:]])
    return times, ram, y_labels

def sum_times(times):
    summed_times = [times[0]]
    i = 1
    for i in range(1,len(times)):
        summed_times.append([x+y for (x,y) in zip(summed_times[i-1], times[i])])
        i+=1
    return summed_times

def main():
    fpath1 = sys.argv[1]
    fpath2 = sys.argv[2]
    out_png = sys.argv[3]

    x_labels = ["MALVA", "VarGeno", "discoSnp++", "BCFtools", "GATK"]
    times1, ram1, y_labels = parse_csv(fpath1)
    times2, ram2, _ = parse_csv(fpath2)
    summed_times1 = sum_times(times1)
    summed_times2 = sum_times(times2)

    index = np.array([0.03, 0.06, 0.09, 0.12, 0.15])
    bar_width = 0.01
    opacity = 1.0
    colors1 = ['blue', 'salmon', 'gold', 'grey']
    colors2 = ['slateblue', 'burlywood', 'khaki', 'silver']

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    x = np.linspace(0,0.18)
    y = x*0
    ax1.plot(x,y,color='black')

    for i in range(0,len(y_labels))[::-1]:
        d1 = summed_times1[i]
        d2 = summed_times2[i]
        ax1.bar(index-bar_width/2, d1, bar_width,
                alpha=opacity,
                color=colors1[i],
                edgecolor='black',
                linewidth=0.5,
                label=y_labels[i] + " (FG)")
        ax2.bar(index+bar_width/2, [-x for x in d2], bar_width,
                alpha=opacity,
                color=colors2[i],
                edgecolor='black',
                linewidth=0.5,
                label=y_labels[i] + " (HG)")

    ax1.bar(index-bar_width/2, [-r for r in ram1], bar_width,
            alpha=opacity,
            color=colors1[-1],
            edgecolor='black',
            linewidth=0.5,
            label="RAM (FG)")
    ax2.bar(index+bar_width/2, ram2, bar_width,
            alpha=opacity,
            color=colors2[-1],
            edgecolor='black',
            linewidth=0.5,
            label="RAM (HG)")

    ax1.set_xlim(-0.02,0.20)
    ax1.set_xticks(index)
    ax1.set_xticklabels(x_labels)
    # ax1.set_xlabel("Tool")
    ax1.set_ylabel("Time (hours)", y = 0.72)
    ax1.set_ylim(-60,80)
    ax1.set_yticks([0,20,40,60,80])
    ax1.legend(loc=2, bbox_to_anchor=(0.003, 0.995))

    ax2.set_ylabel("RAM (GB)", labelpad = 19, rotation=270, y = 0.22)
    ax2.set_ylim(60,-80)
    ax2.set_yticks([0,20,40,60])
    ax2.legend(loc=2, bbox_to_anchor=(0.003, 0.79))

    ax1.grid(color='gray', alpha=0.75, linestyle='dashed', axis='y')
    ax2.grid(color='gray', alpha=0.75, linestyle='dashed', axis='y')

    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)

    plt.subplots_adjust(top=0.97, bottom=0.05, right=0.95, left=0.06)
    DPI = fig.get_dpi()
    w,h = 1366.0, 768.0
    fig.set_size_inches(w/float(DPI),h/float(DPI))
    fig.savefig(out_png, dpi=DPI)
    #plt.show()

if __name__ == '__main__':
    main()
