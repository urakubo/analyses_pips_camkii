from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pprint

black_red = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        }

white_red = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        }
        
black_green = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        ],
        "blue": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        }

white_green = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        }

black_blue = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        ],
        }

white_blue = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        }

black_magenta = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.8906, 0.8906),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 0.0, 0.0),
        (1.0, 0.5, 0.5),
        ],
        }

white_magenta = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 0.8906, 0.8906),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 0.5, 0.5),
        ],
        }

black_cyan = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        ],
        "blue": [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        ],
        }

white_cyan = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        }


white_green_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 3/255, 3/255),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 175/255, 175/255),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 122/255, 122/255),
        ],
        }

white_red_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 75/255, 75/255),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        }
        

white_purple_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 153/255, 153/255),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 153/255, 153/255),
        ],
        }

white_skyblue_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (1.0, 77/255, 77/255),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (1.0, 196/255, 196/255),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        }


white_green_binary_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 3/255),
        (1.0, 3/255, 3/255),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 175/255),
        (1.0, 175/255, 175/255),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 122/255),
        (1.0, 122/255, 122/255),
        ],
        }

white_red_binary_universal = {
    	"red": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        ],
        "green": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 75/255),
        (1.0, 75/255, 75/255),
        ],
        "blue": [
        (0.0, 1.0, 1.0),
        (0.5, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        }
        


c0 = [216,242,85] # Light green
c1 = [255,241,0] # Yellow
c2 = [246,170,0] # Orange
c3 = [255,202,191] # Light pink

gy = 230
c0 = [gy ,gy, gy] # Gray
c1 = [255,202,191] # Ligh pink
c2 = [255,255,128] # Cream
c3 = [gy, gy, gy] # Gray

gy = 230
c0 = [gy ,gy, gy] # Gray
c1 = [255,128,130] # Pink
c2 = [255,202,191] # Ligh pink
c3 = [gy, gy, gy] # Gray


phase_diagram1 = {
    	"red": [
        [0.00, 0.00 , c0[0]],
        [0.25, c0[0], c1[0]],
        [0.50, c1[0], c2[0]],
        [0.75, c2[0], c3[0]],
        [1.00, c3[0], c3[0]],
        ],
        "green": [
        [0.00, 0.00 , c0[1]],
        [0.25, c0[1], c1[1]],
        [0.50, c1[1], c2[1]],
        [0.75, c2[1], c3[1]],
        [1.00, c3[1], c3[1]],
        ],
        "blue": [
        [0.00, 0.00 , c0[2]],
        [0.25, c0[2], c1[2]],
        [0.50, c1[2], c2[2]],
        [0.75, c2[2], c3[2]],
        [1.00, c3[2], c3[2]],
        ],
        }

phase_diagram1 = {k: [[c[0], c[1]/255, c[2] / 255] for c in phase_diagram1[k] ] for k in phase_diagram1}


cmap_phase_diagram1 = LinearSegmentedColormap('phase_diagram1', phase_diagram1, N=256)
cmap_phase_diagram1.set_bad(color='k')


cmap_black_red = LinearSegmentedColormap("black_red", black_red, N=256)
cmap_white_red = LinearSegmentedColormap("white_red", white_red, N=256)
cmap_black_green = LinearSegmentedColormap("black_green", black_green, N=256)
cmap_white_green = LinearSegmentedColormap("white_green", white_green, N=256)
cmap_black_blue = LinearSegmentedColormap("black_blue", black_blue, N=256)
cmap_white_blue = LinearSegmentedColormap("white_blue", white_blue, N=256)
cmap_black_magenta = LinearSegmentedColormap("black_magenta", black_magenta, N=256)
cmap_white_magenta = LinearSegmentedColormap("white_magenta", white_magenta, N=256)
cmap_black_cyan    = LinearSegmentedColormap("black_cyan", black_cyan, N=256)
cmap_white_cyan    = LinearSegmentedColormap("white_cyan", white_cyan, N=256)

cmap_white_green_universal = LinearSegmentedColormap("white_green_universal", white_green_universal, N=256)
cmap_white_red_universal = LinearSegmentedColormap("white_red_universal", white_red_universal, N=256)
cmap_white_purple_universal = LinearSegmentedColormap("white_purple_universal", white_purple_universal, N=256)
cmap_white_skyblue_universal = LinearSegmentedColormap("white_skyblue_universal", white_skyblue_universal, N=256)

cmap_white_green_binary_universal = LinearSegmentedColormap("white_green_binary_universal", white_green_binary_universal, N=256)
cmap_white_red_binary_universal = LinearSegmentedColormap("white_red_binary_universal", white_red_binary_universal, N=256)


cmap = {'CaMKII': cmap_black_green,
		'STG'   : cmap_black_red,
		'GluN2B': cmap_black_magenta,
		'PSD95' : cmap_black_cyan}

cmap_universal_ratio = {
		'CaMKII': (3/255,175/255,122/255),\
		'STG'   : (255/255,75/255,0),\
		'GluN2B': (153/255,0,153/255),\
		'PSD95' : (77/255,196/255,255/255),\
		'All'	: (0.7,0.7,0.7)}

cmap_universal_ratio_light = \
		{k: ((v[0]+2)/3, (v[1]+2)/3, (v[2]+2)/3) for k, v in cmap_universal_ratio.items()}
		
cmap_universal_uint = {
		'CaMKII': (3,175,122),\
		'STG'   : (255,75,0),\
		'GluN2B': (153,0,153),\
		'PSD95' : (77,196,255),\
		'All'	: (177,177,177)}

light_green_universal_uint = (119, 217, 168)
yellow_universal_uint      = (255, 241, 0)

# https://pystyle.info/matplotlib-master-of-colormap/

