from matplotlib import cm



#resolution
dpi=300
alpha_val = 1

#coloring
qualitative_tab_list = ["brown", "orange", "olive", "green", "blue", "cyan", "purple", "pink", "red"]
qualitative_90s_list= ["navy", "blue", "royalblue","darkgreen", "forestgreen", "firebrick", "salmon"]
gradient_kays_list = ['gold', 'orange', 'darkorange', 'tomato', 'orangered', 'red', 'crimson']
gradient_blue_list = ["deepskyblue", "skyblue", "steelblue", "cornflowerblue", "royalblue", "mediumblue", "midgnightblue"]
gradient_green_list = ["chartreuse", "lawngreen", "limegreen", "forestgreen", "seagreen", "green", "darkgreen"]

#maps:
gradient_green_map = cm.get_cmap("Greens")
gradient_blueGreenYellow_map = cm.get_cmap("viridis")
gradient_kays_map = cm.get_cmap("inferno")
qualitative_tab_map = cm.get_cmap("tab20")
qualitative_dark_map = cm.get_cmap("Dark2")



#### ACTIVE STYLE:  ###
active_gradient_map = gradient_green_map
active_qualitative_map = qualitative_tab_map
active_gradient_list = gradient_green_list
active_qualitative_list = qualitative_tab_list
