import matplotlib.pyplot as plt
import pandas as pd
from . import plot_config as pc
from collections import OrderedDict as od

# in_dict = {
# 'group1': [1, 2, 3, 4],
# 'group2': [11, 12, 13, 14],
# 'group3': [21, 22, 23, 24]
# }
# every group is a strip

# group_name_list: from bottom to top
def per_area_plot(out_fig='', in_dict={}, group_name_list=[], subgroup=[], title='', xticks=[], xlabel='', ylabel='Percent (%)'):
    plt.figure()
    # make dataframe from dict
    od_dict = od([(i, in_dict[i]) for i in group_name_list])
    df = pd.DataFrame.from_dict(od_dict)
    # print df

    # transform to percentage dataframe
    perc_df = df.divide(df.sum(axis=1), axis=0) * 100

    # make the plot
    color = pc.get_color_for_group(subgroup)
    perc_df.plot.area(ylim=(0, 100), color=color)

    plt.legend(labelspacing=-2.5, bbox_to_anchor=((1, 0.5, 0.5, 1)), frameon=False, loc=3, ncol=1, mode='expand',
               prop=pc.legend_font)
    plt.margins(0, 0), plt.title(title, **(pc.title_font))
    plt.xlabel(xlabel, **(pc.label_font)), plt.ylabel(ylabel, **(pc.label_font))
    plt.xticks(range(0, len(xticks)), xticks)

    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    # plt.show()


if __name__ == '__main__':
    in_dict = {
        'group1': [11, 12, 13, 14],
        'group2': [11, 12, 13, 14],
        'group3': [21, 22, 23, 24],
        'group4': [11, 12, 13, 14],
        'group5': [21, 22, 23, 24],
    }

    subgroup = [1, 3, 1]
    # color = ['black', 'grey', 'green']
    title = 'Stacked plot'
    xlabel = 'Count >='
    xticks = '12345'
    group_name = ['group1', 'group2', 'group3', 'group4', 'group5']
    out_fig = 'plot.png'
    per_area_plot(title=title, subgroup=subgroup, group_name_list=group_name, out_fig=out_fig, in_dict=in_dict, xlabel=xlabel, xticks=xticks)
