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

black_light_cyan = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.5, 0.5),
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

black_dark_cyan = {
    	"red": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        ],
        "green": [
        (0.0, 0.0, 0.0),
        (1.0, 0.7, 0.7),
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
c0 = [gy ,gy, gy]  # Gray
c1 = [255,128,130] # Pink
c2 = [255,255,128] # Cream
c3 = [0, 90, 255]  # Blue
c3 = [77, 196, 255]  # Light blue



#
def make_colormap_phase_diagram(c):
	ratios = np.linspace(0, 1, len(c)-1, endpoint=True)
	phase_diagram = {
    	"red"  : [[ r, c[i][0]/255, c[i+1][0]/255] for i, r in enumerate(ratios) ],
        "green": [[ r, c[i][1]/255, c[i+1][1]/255] for i, r in enumerate(ratios) ],
        "blue" : [[ r, c[i][2]/255, c[i+1][2]/255] for i, r in enumerate(ratios) ]
        }
	phase_diagram = LinearSegmentedColormap('phase_diagram', phase_diagram, N=256)
	phase_diagram.set_bad(color='k')
	return phase_diagram
#
#
gy   = 230
c    = {}
c[0] = [0.0 , 0.0, 0.0]
c[1] = [gy ,gy, gy] # Gray
c[2] = [255,128,130] # Pink
c[3] = [255,255,255] # Ligh pink
# c[3] = [255,202,191] # Ligh pink
c[4] = [77  ,196, 255]  # Light blue
c[5] = [gy, gy, gy] # Gray
cmap_phase_diagram1 = make_colormap_phase_diagram(c)
#
#
c    = {}
c[0] = [255 ,128, 130] # Pink
c[1] = [255 ,128, 130] # Pink
c[2] = [255 ,255, 255] # White
c[3] = [77  ,196, 255]  # Light blue
c[4] = [gy  , gy,  gy] # Gray
c[5] = [gy  , gy,  gy] # Gray
cmap_phase_diagram2 = make_colormap_phase_diagram(c)
#
#
c    = {}
c[0] = [255 ,255, 255] # White
c[1] = [255 ,255, 255] # White
c[2] = [gy  , gy,  gy] # Gray
c[3] = [77+89  ,196+30, 255]  # Light blue
c[4] = [77+89  ,196+30, 255]  # Light blue
cmap_phase_diagram3 = make_colormap_phase_diagram(c)



c    = {}
c[0] = [0.0 , 0.0, 0.0]
c[1] = [gy ,gy, gy] # Gray
c[2] = [255,128,130] # Pink
#c[3] = [255,255,255] # Ligh pink
# c[3] = [255,202,191] # Ligh pink
c[3] = [77  ,196, 255]  # Light blue
c[4] = [77  ,196, 255]  # Light blue
cmap_phase_diagram4 = make_colormap_phase_diagram(c)

c    = {}
c[0] = [255 ,255, 255] # White
c[1] = [255 ,255, 255] # White
c[2] = [255 ,255, 255] # White
c[3] = [255 ,255, 255] # White
c[4] = [255 ,255, 255] # White
cmap_phase_diagram5 = make_colormap_phase_diagram(c)

c    = {}
c[0] = [255,191,192] # Light pink
c[1] = [255,191,192] # Light pink
c[2] = [255,191,192] # Light pink
c[3] = [255,191,192] # Light pink
c[4] = [255,191,192] # Light pink
cmap_phase_diagram6 = make_colormap_phase_diagram(c)


c    = {}
c[0] = [0.0 , 0.0, 0.0]
c[1] = [gy ,gy, gy] # Gray
c[2] = [77  ,196, 255]  # Light blue
c[3] = [77  ,196, 255]  # Light blue
cmap_phase_diagram7 = make_colormap_phase_diagram(c)


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
cmap_black_light_cyan    = LinearSegmentedColormap("black_light_cyan", black_light_cyan, N=256)
cmap_black_dark_cyan    = LinearSegmentedColormap("black_dark_cyan", black_dark_cyan, N=256)

cmap_white_green_universal = LinearSegmentedColormap("white_green_universal", white_green_universal, N=256)
cmap_white_red_universal = LinearSegmentedColormap("white_red_universal", white_red_universal, N=256)
cmap_white_purple_universal = LinearSegmentedColormap("white_purple_universal", white_purple_universal, N=256)
cmap_white_skyblue_universal = LinearSegmentedColormap("white_skyblue_universal", white_skyblue_universal, N=256)

cmap_white_green_binary_universal = LinearSegmentedColormap("white_green_binary_universal", white_green_binary_universal, N=256)
cmap_white_red_binary_universal = LinearSegmentedColormap("white_red_binary_universal", white_red_binary_universal, N=256)


cmap = {'CaMKII': cmap_black_green,
		'STG'   : cmap_black_red,
		'GluN2B': cmap_black_magenta,
		'PSD95' : cmap_black_cyan,
		'All PSD95' : cmap_black_cyan,
		'Shared PSD95' : cmap_black_dark_cyan,
		'Unshared PSD95' : cmap_black_light_cyan,
		'PSD95 shared only by GluN2B' : cmap_black_dark_cyan}

cmap_universal_ratio = {
		'CaMKII': (3/255,175/255,122/255),\
		'STG'   : (255/255,75/255,0),\
		'GluN2B': (153/255,0,153/255),\
		'PSD95' : (77/255,196/255,255/255),\
		'Shared PSD95' : (0/255,90/255,255/255),
		'Unshared PSD95' : (191/255,228/255,255/255),\
		'CaMKII hub'	: (0.3,0.3,0.3),\
		'All'	: (0.3,0.3,0.3),\
		'PIPS': (60/255,60/255,60/255)}

cmap_universal_ratio_light = \
		{k: ((v[0]+2)/3, (v[1]+2)/3, (v[2]+2)/3) for k, v in cmap_universal_ratio.items()}
		
cmap_universal_uint = {
		'CaMKII': (3,175,122),\
		'STG'   : (255,75,0),\
		'GluN2B': (153,0,153),\
		'PSD95' : (77,196,255),\
		'Shared PSD95' : (0,90,255),
		'Unshared PSD95' : (191,228,255),\
		'All'	: (177,177,177)}

light_green_universal_uint = (119, 217, 168)
light_green_universal_ratio = (119/255, 217/255, 168/255)

ellow_universal_uint      = (255, 241, 0)

green_universal_uint = (3,175,122)

red_uint     = (255,75,0)
green_uint   = (3,175,122)
blue_uint    = (0,90,255)
pink_uint    = (255,128,130)
orange_uint  = (246,170,0)
skyblue_uint = (77,196,255)
purple_uint  = (153,0,153)
blown_uint   = (128,64,0)
yellow_uint  = (255,241,0)
cols_list_uint  = [\
	red_uint, green_uint,  blue_uint, pink_uint, orange_uint, skyblue_uint, purple_uint, blown_uint, yellow_uint]
	
cols_list_ratio = [[cols[0]/255, cols[1]/255, cols[2]/255] for cols in cols_list_uint]

#col_list_dendro = 
# https://pystyle.info/matplotlib-master-of-colormap/


