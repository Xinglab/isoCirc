import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from . import plot_config as pc


def stackedBarPlot_core(out_fig='', in_dict={}, subgroup=[], title='', ystep=200, xticks=[], xlabel='', ylabel=''):
    plt.figure()
    # make dataframe
    df = pd.DataFrame.from_dict(in_dict)
    y_max = max(in_dict.values()[0])
    # make line plot
    color = pc.get_color_for_group(subgroup)
    df.plot(kind='bar', stacked=True)  # , color=color) # use default colors

    plt.margins(0, 0), plt.title(title, **(pc.title_font))
    plt.xlabel(xlabel, **(pc.label_font)), plt.ylabel(ylabel, **(pc.label_font))

    # yticks
    y_max = max(in_dict.values()[0])
    pc.set_y_ticks(y_max, plt)

    # xticks
    plt.xticks(range(0, len(xticks)), xticks)
    plt.tick_params(axis='x', labelrotation=45)  # labelsize=8,
    plt.legend(labelspacing=-2.5, bbox_to_anchor=((1, 0.5, 0.5, 1)), frameon=False, loc=3, ncol=1, mode='expand',
               prop=pc.legend_font)

    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    # plt.show()


if __name__ == '__main__':
    in_dict = {
        'group1': [5, 10, 20],
        'group2': [10, 5, 13],
        'group3': [15, 6, 4]
    }
    subgroup = [1, 1, 1]
    out_fig = './box.png'
    xticks = [1, 2, 3]
    stackedBarPlot_core(out_fig=out_fig, in_dict=in_dict, subgroup=subgroup, xticks=xticks)
