import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib import rc
rc('text')
font = {'family' : 'normal',
        'size'   : 17}
matplotlib.rc('font', **font)

def main():
    fpath = sys.argv[1]
    out_png = sys.argv[2]

    Rsnp = {}
    Psnp = {}
    Rindel = {}
    Pindel = {}
    for line in open(fpath, 'r'):
        if line[0] == 'T':
            continue
        tool,chrom,rsnp,psnp,rindel,pindel = line.strip("\n").split(",")
        rsnp,psnp,rindel,pindel = round(float(rsnp),3),round(float(psnp),3),round(float(rindel),3),round(float(pindel),3)
        Rsnp[tool] = Rsnp[tool] + [rsnp] if tool in Rsnp else [rsnp]
        Psnp[tool] = Psnp[tool] + [psnp] if tool in Psnp else [psnp]
        Rindel[tool] = Rindel[tool] + [rindel] if tool in Rindel else [rindel]
        Pindel[tool] = Pindel[tool] + [pindel] if tool in Pindel else [pindel]

    tools = ["malva", "vargeno", "discosnp", "bcftools", "gatk"] #sorted(Pindel.keys())
    labels = ["MALVA", "VarGeno", "DiscoSnp", "BCFtools", "GATK"]

    vg = "vargeno"
    if vg not in Rsnp:
        Rsnp[vg] = [0]
        Psnp[vg] = [0]
        Rindel[vg] = [0]
        Pindel[vg] = [0]
    
    psnp_color = 'darkgreen'
    rsnp_color = 'lightgreen'
    pindel_color = 'crimson'
    rindel_color = 'coral'
    alpha = 1
    edge_color = 'black'
    edge_size = 0.5
    d = 3
    base_positions_1 = [1+d*x for x in range(0,5)]
    base_positions_2 = [x+d*6 for x in base_positions_1]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    Ps_vp = ax1.violinplot([Psnp[x] for x in tools], positions=base_positions_1, showextrema=False)
    for pc in Ps_vp['bodies']:
        pc.set_facecolor(psnp_color)
        pc.set_edgecolor(edge_color)
        pc.set_linewidths(edge_size)
        pc.set_alpha(alpha)
    
    Rs_vp = ax1.violinplot([Rsnp[x] for x in tools], positions=[x+0.5 for x in base_positions_1], showextrema=False)
    for pc in Rs_vp['bodies']:
        pc.set_facecolor(rsnp_color)
        pc.set_edgecolor(edge_color)
        pc.set_linewidths(edge_size)
        pc.set_alpha(alpha)

    Pi_vp = ax1.violinplot([Pindel[x] for x in tools], positions=base_positions_2, showextrema=False)
    for pc in Pi_vp['bodies']:
        pc.set_facecolor(pindel_color)
        pc.set_edgecolor(edge_color)
        pc.set_linewidths(edge_size)
        pc.set_alpha(alpha)
    
    Rs_vp = ax1.violinplot([Rindel[x] for x in tools], positions=[x+0.5 for x in base_positions_2], showextrema=False)
    for pc in Rs_vp['bodies']:
        pc.set_facecolor(rindel_color)
        pc.set_edgecolor(edge_color)
        pc.set_linewidths(edge_size)
        pc.set_alpha(alpha)

    # plt.title("Results on full dataset")
    ax1.set_xlim([0,base_positions_2[-1] + 1.5])
    ax1.set_xticks([(base_positions_1[0]+base_positions_1[-1])/2+0.25, (base_positions_2[0]+base_positions_2[-1])/2+0.25])
    ax1.set_xticklabels(["SNP", "INDEL"])
    ax1.set_ylim([0,1])
    ax1.set_yticks([0,0.25,0.5,0.75,1])

    ax2.set_xlim([0,base_positions_2[-1] + 1.5])
    ax2.set_xticks([x+0.25 for x in base_positions_1 + base_positions_2])
    ax2.set_xticklabels(labels + labels)

    

    patch1 = mpatches.Patch(facecolor=psnp_color, edgecolor=edge_color, linewidth=edge_size, label='SNP (precision)')
    patch2 = mpatches.Patch(facecolor=rsnp_color, edgecolor=edge_color, linewidth=edge_size, label='SNP (recall)')
    patch3 = mpatches.Patch(facecolor=pindel_color, edgecolor=edge_color, linewidth=edge_size, label='Indel (precision)')
    patch4 = mpatches.Patch(facecolor=rindel_color, edgecolor=edge_color, linewidth=edge_size, label='Indel (recall)')
    ax1.legend(handles=[patch1, patch2], loc=3)
    ax2.legend(handles=[patch3, patch4], loc=4)

    ax1.grid(color='gray', alpha=0.75, linestyle='dashed', axis='y')
    ax2.grid(color='gray', alpha=0.75, linestyle='dashed', axis='x')

    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)

    plt.subplots_adjust(top=0.95, bottom=0.03, right=0.99, left=0.05)
    DPI = fig.get_dpi()
    fig.set_size_inches(1280.0/float(DPI),1024.0/float(DPI))
    fig.savefig(out_png, dpi=DPI)
    # plt.show()

if __name__ == "__main__":
    main()
