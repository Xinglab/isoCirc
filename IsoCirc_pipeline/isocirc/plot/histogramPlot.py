import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from . import plot_config as pc


def histogramPlot_core(out_fig='', bins=[], in_dict={}, title='', xticks=[], xlabel='', ylabel=''):
    plt.figure()
    # make dataframe
    df = pd.DataFrame.from_dict(in_dict)
    # make line plot
    plt.hist(df.values, bins=bins, edgecolor='k', color='g')  # use default color
    # plt.show()

    plt.margins(0, 0), plt.title(title, **(pc.title_font))
    plt.xlabel(xlabel, **(pc.label_font)), plt.ylabel(ylabel, **(pc.label_font))

    bin_cnt = []
    for i in range(0, len(bins) - 1):
        bin_cnt.append(0)
        for x in sum(in_dict.values(), []):
            if bins[i] <= x < bins[i + 1]:
                bin_cnt[-1] += 1
    # yticks
    y_max = max(bin_cnt)
    pc.set_y_ticks(y_max, plt)

    plt.xticks(bins)
    plt.tick_params(axis='x', labelsize=8, labelrotation=45)
    # plt.legend(labelspacing=-2.5, bbox_to_anchor=((1, 0.5, 0.5, 1)), frameon=False, loc=3, ncol=1, mode='expand', prop=pc.legend_font)

    plt.savefig(out_fig, dpi=300, bbox_inches="tight")


if __name__ == '__main__':
    in_dict = {
        'group1': [5, 6, 5, 6, 6, 7, 10, 20, 30, 1, 2],
    }
    out_fig = './box.png'
    bins = [0, 20, 40]
    xticks = bins
    stackedBarPlot(out_fig=out_fig, bins=bins, in_dict=in_dict, xticks=xticks, title='barPlot',
                   xlabel='count', ylabel='isoformCnt')
