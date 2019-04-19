# color
red_color = ['red', 'pink', 'tomato', 'firebrick', 'salmon', 'darkred', ]
blue_color = ['steelblue', 'blue', 'darkblue', 'cyan', 'skyblue', 'dodgerblue']
green_color = ['green', 'lightgreen', 'greenyellow', 'seagreen',  'limegreen', 'mediumaquamarine']
yellow_color = ['orange', 'yellow', 'gold', 'chocolate', 'goldenrod', 'peachpuff']
black_color = ['black', 'gray', 'silver', 'white']
whole_color = [blue_color, yellow_color, red_color, green_color, black_color]

def get_color_for_group(group=[]):
    color = []
    for i, n in enumerate(group):
        color.extend(whole_color[i][:n])
    return color

step_list = []
for i in range(6):
    step_list.extend([0.1 * pow(10, i), 0.2 * pow(10, i), 0.5 * pow(10, i)])
ticks_num = [5, 6, 7, 8, 9]

def set_y_ticks(max_value, plt):
    y_step, y_num = -1, -1
    for step in step_list:
        for num in ticks_num:
            if step * num >= max_value:
                y_step = int(step) if step > 1 else step
                y_num = num
                break
        if y_step != -1:
            break
    if y_step == -1:
        return [], [], -1

    sub_yticks_idx = range(0, y_num + 1)
    sub_yticks = [i * y_step for i in range(0, y_num + 1)]
    # set y ticks
    plt.ylim((0, y_step * y_num))
    plt.yticks(sub_yticks, sub_yticks)
    return


# font
# {'fontname': 'sans-serif', 'fontsize':10, 'fontweight':'bold', 'fontstyle':'italic'}
title_font = {'fontname': 'sans', 'fontsize': 'x-large', 'fontweight': 'bold', 'fontstyle': 'normal'}
# label_font = {'fontname': 'sans-serif', 'fontsize': 'large', 'fontweight': 'normal', 'fontstyle': 'normal'}
# legend_font = {'fontname': 'sans-serif', 'fontsize': 'medium', 'fontweight': 'normal', 'fontstyle': 'normal'}

# title_font = {'name': 'sans-serif', 'size': 'x-large', 'weight': 'bold', 'style': 'normal'}
label_font = {'family': 'sans-serif', 'size': 'large', 'weight': 'bold', 'style': 'normal'}
legend_font = {'family': 'sans-serif', 'size': 'medium', 'weight': 'bold', 'style': 'normal'}
