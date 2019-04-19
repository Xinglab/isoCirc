import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from . import plot_config as pc


def dis_to_known_plot(out_fig='', in_dict={}, subgroup=[], title='', xticks=[], xlabel='', ylabel=''):
    plt.figure()
    # make dataframe
    df = pd.DataFrame.from_dict(in_dict)
    # make line plot
    color = pc.get_color_for_group(subgroup)
    df.plot.line(color=color, linewidth=2)

    plt.margins(0, 0), plt.title(title, **(pc.title_font))
    plt.xlabel(xlabel, **(pc.label_font)), plt.ylabel(ylabel, **(pc.label_font))

    # yticks
    y_max = max(sum(in_dict.values(), []))
    pc.set_y_ticks(y_max, plt)
    # xticks
    sub_xticks_idx, sub_xticks = [], []
    xstep = 5
    for i in range(len(xticks) / xstep + 1):
        if i * xstep < len(xticks):
            sub_xticks_idx.append(i*xstep)
            sub_xticks.append(xticks[i*xstep])

    # plt.xticks(range(0, len(xticks)), xticks)
    plt.xticks(sub_xticks_idx, sub_xticks)
    plt.legend(labelspacing=-2.5, bbox_to_anchor=((1, 0.5, 0.5, 1)), frameon=False, loc=3, ncol=1, mode='expand',
               prop=pc.legend_font)

    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    # plt.show()


if __name__ == '__main__':
    in_dict = {
        'group1': [11, 12, 13, 14, 15],
        'group2': [12, 14, 12, 11, 12],
        'group3': [21, 22, 23, 24, 20],
        'group4': [22, 20, 24, 23, 25],
    }
    subgroup = [2, 2]
    out_fig = './dis.png'
    xticks = [-2, -1, 0, 1, 2]
    disToKnownPlot_core(out_fig=out_fig, in_dict=in_dict, subgroup=subgroup, xticks=xticks)
