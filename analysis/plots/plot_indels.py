import sys

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from matplotlib import rc
rc('text')
font = {'family' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

def main():
    # gt = sys.argv[1]
    lmin = -100 # int(sys.argv[2])
    lmax = 200 # int(sys.argv[3])
    fpaths = sys.argv[1:-1]
    out_pdf = sys.argv[-1]

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212) # ax1.twinx()

    tools = ['MALVA', 'GATK', 'BCFtools', 'discoSnp++']
    colors = ['red', 'green', 'blue', 'orange']
    i = 0
    for fpath in fpaths:
        tps = {}
        fps = {}
        tots = {}
        for line in open(fpath):
            if line[0:4] == 'chr,':
                continue
            chrom, _gt, l, tp, fp, tot = line.strip('\n').split(',')
            if tot == "0": # These are FPs, we don't need them
                continue
            # if _gt != gt:
            #     continue
            l = int(l)
            if lmin <= l <= lmax:
                tps[l] = tps[l] + int(tp) if l in tps else int(tp)
                fps[l] = fps[l] + int(fp) if l in fps else int(fp)
                tots[l] = tots[l] + int(tot) if l in tots else int(tot)

        tots_mod = {l:tots[l]+1 for l in tots}
        ax1.scatter(sorted(tps.keys()), [tps[l]/tots[l] for l in sorted(tps.keys())], color=colors[i], label=tools[i], linewidths=0.0001, alpha=0.75, s=23)
        if i == 0:
            ax2.bar(sorted(tps.keys()), [np.log(tots_mod[l]) for l in sorted(tps.keys())], color="grey")
        i += 1
    # plt.xticks(np.arange(min(Xs), max(Xs)+1, 25))

    ax1.legend(loc=4, bbox_to_anchor=(1, -0.17), ncol=4)
    #ax1.set_title("Recall on {} indels".format(gt))
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel("Recall")
    ax2.set_xlabel("Indel length (#bp)")
    ax2.set_ylabel("#indels (log scale)")
    ax2.set_ylim(0,12)
    # xlabel('Item (s)')
    # ylabel('Value')
    # title('Python Line Chart: Plotting numbers')
    # grid(True)
    plt.subplots_adjust(top=0.99, bottom=0.09, right=0.99, left=0.07)
    DPI = fig.get_dpi()
    fig.set_size_inches(1366.0/float(DPI),768.0/float(DPI))
    fig.savefig(out_pdf, dpi=DPI)
    # plt.show()

if __name__ == '__main__':
    main()
