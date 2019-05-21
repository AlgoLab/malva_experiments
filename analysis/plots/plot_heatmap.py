import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns; sns.set()

# from matplotlib import rc
# rc('text', usetex=True)

# font = {'family' : 'normal',
#         'size'   : 15}
# matplotlib.rc('font', **font)

sns.set(font_scale=1.3)

category_labels = ['HomoRef', 'HetRef', 'HomoAlt', 'HetAlt']

def norm_table_col(table):
    norm_table = []
    colsum = [sum(x) for x in zip(*table)]
    colsum[-1] = 0
    for row in table:
        norm_table.append([round((x*1.0)/y, 3) if y != 0 else 0 for x,y in zip(row, colsum)])
    return norm_table

def norm_table_row(table):
    norm_table = []
    rowsum = [sum(row) for row in table]
    i = 0
    for row in table:
        s = rowsum[i]
        norm_table.append([round((x*1.0)/s, 3) if s != 0 else 0 for x in row])
        i+=1
    return norm_table

def draw_table(table, labels, out_prefix):
    ax = sns.heatmap(table,
                     vmin = 0,
                     vmax = 0, # no norm
                     #vmax = 1,
                     annot=True,
                     fmt='d', # no norm
                     #fmt='.3f',
                     linecolor = 'black', # no norm
                     linewidths = 1,
                     cbar = False, # no norm
                     cmap="Greys",
                     square = True,
                     xticklabels=labels + ["Uncalled"],
                     yticklabels=labels)
    ax.invert_yaxis()
    #plt.title(r"\textbf{MALVA - SNPs Analysis on Half-DS}")
    plt.xlabel(r"\textbf{Given GT}")
    plt.ylabel(r"\textbf{Real GT}")
    plt.yticks(np.arange(len(labels))+0.5,labels, rotation=90, va="center")
    # plt.show()
    plt.savefig(out_prefix + ".table.png")

def draw_heatmap(table, labels, out_prefix):
    ax = sns.heatmap(table,
                     vmin = 0,
                     vmax = 1,
                     annot=True,
                     fmt='.3f',
                     linewidths = 1,
                     cmap="Greys",
                     square = True,
                     xticklabels=labels + ["Uncalled"],
                     yticklabels=labels)
    ax.invert_yaxis()
    #plt.title(r"\textbf{MALVA - SNPs Analysis on Half-DS}")
    plt.xlabel(r"\textbf{Given GT}")
    plt.ylabel(r"\textbf{Real GT}")
    plt.yticks(np.arange(len(labels))+0.5,labels, rotation=90, va="center")
    #plt.show()
    plt.savefig(out_prefix + ".heatmap.png")

def get_gt_category(gt):
    # category_labels = ['HomoRef', 'HetRef', 'HomoAlt', 'HetAlt']
    (x,y) = (int(x) for x in gt.split('/'))
    if x == y and x == 0:
        return 0 # HomoRef
    elif x == y:
        return 2 # HomoAlt
    elif x != y and (x == 0 or y == 0):
        return 1 # HetRef
    else:
        return 3 # HetAlt

def parse_csv(csv_path):
    table = []
    labels = []
    header_flag = True
    for line in open(csv_path, 'r'):
        if header_flag:
            labels = line.strip('\n').split(',')[1:-1] #-1 for Uncalled
            header_flag = False
            continue
        row = [int(x) for x in line.strip('\n').split(',')[1:]]
        table.append(row)
    return table, labels

def build_cumulate_table(table, labels):
    cumulate_table = [[0]*5 for _ in range(4)]
    for i in range(0,len(labels)):
        truth_gt = labels[i]
        truth_cat = get_gt_category(truth_gt)
        for j in range(0,len(labels)):
            given_gt = labels[j]
            given_cat = get_gt_category(given_gt)
            value = int(table[i][j])
            cumulate_table[truth_cat][given_cat] += value
        j = len(labels)
        value = int(table[i][j])
        cumulate_table[truth_cat][4] += value
    return cumulate_table, category_labels

def main():
    csv_path = sys.argv[1]
    out_prefix = sys.argv[2]

    table, labels = parse_csv(csv_path)

    table, labels = build_cumulate_table(table, labels)
    draw_table(table, labels, out_prefix)

    table = norm_table_row(table)
    draw_heatmap(table, labels, out_prefix)

if __name__ == '__main__':
    main()
