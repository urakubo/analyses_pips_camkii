from matplotlib.colors import LinearSegmentedColormap
import numpy as np

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
		'All'	: (0.0,0.0,0.0)}

cmap_universal_uint = {
		'CaMKII': (3,175,122),\
		'STG'   : (255,75,0),\
		'GluN2B': (153,0,153),\
		'PSD95' : (77,196,255),\
		'All'	: (0,0,0)}


# https://pystyle.info/matplotlib-master-of-colormap/

