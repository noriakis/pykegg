import requests
import numpy as np
import re
import os
import matplotlib as mpl
from PIL import Image
import cv2
import pandas as pd
from plotnine import ggplot, geom_point, aes, geom_segment, theme_void, options, geom_rect, geom_text, geom_raster, geom_path


def overlay_opencv_image(nds, 
                         path=None, 
                         pid=None,
                         fill_color="color", 
                         transparent_colors=["#FFFFFF","#BFFFBF"]):
    nds["x0"] = nds["x"] - nds["width"] /2
    nds["y0"] = nds["y"] + nds["height"] /2
        
    if path!=None and os.path.isfile(path):
        im = cv2.imread(path)
    else:
        if pid!=None:
            pid_pref = "".join(re.findall("[a-zA-Z]+", pid))
            im_res = requests.get("https://www.kegg.jp/kegg/pathway/"+pid_pref+"/"+pid+".png")
            nparr = np.frombuffer(im_res.content, np.uint8)
            im = cv2.imdecode(nparr, cv2.IMREAD_COLOR)

    canvas = np.zeros([im.shape[0],im.shape[1],3],dtype=np.uint8)
    canvas.fill(255)
    
    dst = cv2.cvtColor(im, cv2.COLOR_BGR2BGRA)
    for col in transparent_colors:
        h = col[1:7]
        cand_color = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
        mask = np.all(im[:,:,:] == list(cand_color), axis=-1)
        dst[mask,3] = 0
        
    for i in nds.index:
        tmp = nds.iloc[i,:]
        pos = (int(tmp["x0"]),
               int(-1*tmp["y0"]),
               int(tmp["width"]),
               int(tmp["height"]))

        tmp_col = tmp["color"]
        if isinstance(tmp_col, list):
            num_col = len(tmp_col)
            nudge = tmp["width"] / num_col
            for e, one_tmp_col in enumerate(tmp_col):
                new_width = tmp["width"] - (nudge*e)
                canvas = cv2.rectangle(img=canvas,
                           rec=(int(tmp["x0"]),
                            int(-1*tmp["y0"]),
                            int(new_width),int(tmp["height"])),
                           color=hex2rgb(one_tmp_col),
                           thickness=-1)          

        else:
            canvas = cv2.rectangle(img=canvas,
                       rec=pos, color=hex2rgb(tmp_col),
                       thickness=-1)

    image = overlay(canvas, dst)
    return(image)

def overlay(rects, kegg_map):
    rects = cv2.cvtColor(rects, cv2.COLOR_BGR2RGB)
    rects = Image.fromarray(rects).convert("RGBA")
    kegg_map = cv2.cvtColor(kegg_map, cv2.COLOR_BGRA2RGBA)
    kegg_map = Image.fromarray(kegg_map).convert('RGBA')
    result_image = Image.alpha_composite(rects, kegg_map)
    return cv2.cvtColor(np.asarray(result_image), cv2.COLOR_RGBA2BGRA)

def plot_kegg_pathway_plotnine(g, node_x_nudge=20, node_y_nudge=10,
                              split_graphics_name=True, subtype_num=0,
                               label_size=2, show_label="gene",
                              edge_color="subtype", text_label="graphics_name"):
    nddf = g.get_nodes(node_x_nudge=node_x_nudge, node_y_nudge=node_y_nudge)
    eddf = g.get_edges()
    if split_graphics_name:
        nddf["graphics_name"] = nddf.graphics_name.apply(lambda x: x.split(",")[0])

    seg_df = pd.concat([
        nddf.reset_index().set_index("id").loc[eddf.entry1].reset_index().loc[:,["x","y"]],
        nddf.reset_index().set_index("id").loc[eddf.entry2].reset_index().loc[:,["x","y"]]
    ], axis=1)
    seg_df.columns = ["x","y", "xend", "yend"]
    seg_df = pd.concat([seg_df, eddf], axis=1)

    seg_df["subtype"] = seg_df.subtypes.apply(lambda x: x[0][subtype_num] if x!=None else x)
    
    plot = (ggplot() +
     geom_segment(aes(x="x", y="y",xend="xend",yend="yend",color=edge_color), data=seg_df) + 
     geom_rect(aes(xmin="xmin", ymin="ymin", xmax="xmax", ymax="ymax"),
               data=nddf[nddf.original_type=="gene"], fill="white", color="grey") +
     geom_rect(aes(xmin="xmin", ymin="ymin", xmax="xmax", ymax="ymax"),
               data=nddf[nddf.original_type=="compound"], fill="white", color="grey") +
     geom_text(aes(x="x", y="y", label=text_label, filter="original_type!='group'"),
               data=nddf[nddf.original_type==show_label], size=label_size)+theme_void())
    return(plot)

def plot_kegg_global_map_plotnine(g, hide_map=True):
    nds = g.get_nodes()
    coords = g.get_coords()
    if hide_map:
        nds = nds[nds.original_type!="map"]
    plt = (ggplot() +
     geom_segment(aes(x="x",y="y", xend="xend", yend="yend"),
                  color=coords["fgcolor"].tolist(), data=coords) + 
     geom_point(aes(x="x",y="y"), color=nds["fgcolor"].tolist(), data=nds)+
     theme_void())
    return(plt)

def color_grad(minimum=-2,maximum=2,seq=0.1,
              min_col="#ffffff",max_col="#ff0000",
              round_num=2):
    minmax = np.arange(minimum,maximum,seq)
    N = len(minmax)
    conv={}
    for x in range(N):
        l=np.round(minmax[x],round_num)
        c1=np.array(mpl.colors.to_rgb(min_col))
        c2=np.array(mpl.colors.to_rgb(max_col))
        conv[l] = mpl.colors.to_hex((1-x/N)*c1 + x/N*c2)
    return(conv)

def hex2rgb(h):
    return(tuple(int(h.lstrip("#")[i:i+2], 16) for i in (0, 2, 4)))