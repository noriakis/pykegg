import re
import os
import requests
import cv2

import numpy as np
import matplotlib as mpl
import pandas as pd

from PIL import Image
from plotnine import (
    ggplot,
    geom_point,
    aes,
    geom_segment,
    theme_void,
    geom_rect,
    geom_text,
)


def overlay_opencv_image(
        node_df,
        path=None,
        pid=None,
        fill_color="color",
        transparent_colors=None,
):
    """Obtain the raw image of pathway and color the nodes.

    Parameters:
    -----------
    node_df: DataFrame
        node data obtained by `get_nodes()`.
    path: str
        path to the image if already downloaded.
    pid: str
        KEGG pathway identifier.
    fill_color: str
        the column in `node_df` specifying color in HEX.
        If list is given, split the width according to
        the color number.
    transparent_color: list of str
        specify which color to be transparent.
        If `None`, default `["#FFFFFF", "#BFFFBF"]` is used.

    """
    if transparent_colors is None:
        transparent_colors = ["#FFFFFF", "#BFFFBF"]

    node_df["x0"] = node_df["x"] - node_df["width"] / 2
    node_df["y0"] = node_df["y"] + node_df["height"] / 2

    if path is not None and os.path.isfile(path):
        image = cv2.imread(path)
    else:
        if pid is not None:
            pid_pref = "".join(re.findall("[a-zA-Z]+", pid))
            im_res = requests.get(
                "https://www.kegg.jp/kegg/pathway/" + pid_pref + "/" + pid + ".png"
            )
            nparr = np.frombuffer(im_res.content, np.uint8)
            image = cv2.imdecode(nparr, cv2.IMREAD_COLOR)

    canvas = np.zeros([image.shape[0], image.shape[1], 3], dtype=np.uint8)
    canvas.fill(255)

    dst = cv2.cvtColor(image, cv2.COLOR_BGR2BGRA)
    for col in transparent_colors:
        hex_str = col[1:7]
        cand_color = tuple(int(hex_str[i : i + 2], 16) for i in (0, 2, 4))
        mask = np.all(image[:, :, :] == list(cand_color), axis=-1)
        dst[mask, 3] = 0

    for i in node_df.index:
        tmp = node_df.iloc[i, :]
        pos = (
            int(tmp["x0"]),
            int(-1 * tmp["y0"]),
            int(tmp["width"]),
            int(tmp["height"]),
        )

        tmp_col = tmp[fill_color]
        if isinstance(tmp_col, list):
            num_col = len(tmp_col)
            nudge = tmp["width"] / num_col
            for col_num, one_tmp_col in enumerate(tmp_col):
                new_width = tmp["width"] - (nudge * col_num)
                canvas = cv2.rectangle(
                    img=canvas,
                    rec=(
                        int(tmp["x0"]),
                        int(-1 * tmp["y0"]),
                        int(new_width),
                        int(tmp["height"]),
                    ),
                    color=hex2rgb(one_tmp_col),
                    thickness=-1,
                )

        else:
            canvas = cv2.rectangle(
                img=canvas, rec=pos, color=hex2rgb(tmp_col), thickness=-1
            )

    image = overlay(canvas, dst)
    return image


def overlay(rects, kegg_map):
    """overlay two images with transparency

    Parameters:
    -----------
    rects: np.array
    kegg_map: np.array
    """
    rects = cv2.cvtColor(rects, cv2.COLOR_BGR2RGB)
    rects = Image.fromarray(rects).convert("RGBA")
    kegg_map = cv2.cvtColor(kegg_map, cv2.COLOR_BGRA2RGBA)
    kegg_map = Image.fromarray(kegg_map).convert("RGBA")
    result_image = Image.alpha_composite(rects, kegg_map)
    return cv2.cvtColor(np.asarray(result_image), cv2.COLOR_RGBA2BGRA)


def plot_kegg_pathway_plotnine(
        graph,
        node_x_nudge=20,
        node_y_nudge=10,
        split_graphics_name=True,
        subtype_num=0,
        label_size=2,
        show_label="gene",
        edge_color="subtype",
        text_label="graphics_name",
):
    """Plot KEGG pathway using plotnine.
    Parameters:
    -----------
    graph: KGML_graph class
        graph object obtained by `KGML_graph()`.
    node_x_nudge: int
        nudge the x position of nodes.
    node_y_nudge: int
        nudge the y position of nodes.
    split_graphics_name: bool
        if True, split the graphics_name by comma and use the first one.
    subtype_num: int
        the number of subtypes to be used for the edge label.
    label_size: int
        the size of the label.
    show_label: str
        the type of label to be shown, e.g. "gene", "compound".
    edge_color: str
        the column in `edge_df` specifying edge color.
    text_label: str
        the column in `node_df` specifying node label.
    """
    node_df = graph.get_nodes(node_x_nudge=node_x_nudge, node_y_nudge=node_y_nudge)
    edge_df = graph.get_edges()
    if split_graphics_name:
        node_df["graphics_name"] = node_df.graphics_name.apply(lambda x: x.split(",")[0])

    ## Collapse subtypes
    edge_df_col = []
    for i in edge_df.index:
        tmp = edge_df.iloc[i, :]
        for subtype in tmp.subtypes:
            edge_df_col.append([tmp.entry1, tmp.entry2, tmp.type, subtype, tmp.reaction])
    edge_df = pd.DataFrame(edge_df_col)
    edge_df.columns = ["entry1", "entry2", "type", "subtypes", "reaction"]
    seg_df = pd.concat(
        [
            node_df.reset_index()
            .set_index("id")
            .loc[edge_df.entry1]
            .reset_index()
            .loc[:, ["x", "y"]],
            node_df.reset_index()
            .set_index("id")
            .loc[edge_df.entry2]
            .reset_index()
            .loc[:, ["x", "y"]],
        ],
        axis=1,
    )
    seg_df.columns = ["x", "y", "xend", "yend"]
    seg_df = pd.concat([seg_df, edge_df], axis=1)

    ## [TODO] implement multiple edges like `geom_edge_parallel` in ggraph
    ## Currently edges are overlapped.
    seg_df["subtype"] = seg_df.subtypes.apply(
        lambda x: x[subtype_num] if x is not None else x
    )

    plot = (
        ggplot()
        + geom_segment(
            aes(x="x", y="y", xend="xend", yend="yend", color=edge_color), data=seg_df
        )
        + geom_rect(
            aes(xmin="xmin", ymin="ymin", xmax="xmax", ymax="ymax"),
            data=node_df[node_df.original_type == "gene"],
            fill="white",
            color="grey",
        )
        + geom_rect(
            aes(xmin="xmin", ymin="ymin", xmax="xmax", ymax="ymax"),
            data=node_df[node_df.original_type == "compound"],
            fill="white",
            color="grey",
        )
        + geom_text(
            aes(x="x", y="y", label=text_label, filter="original_type!='group'"),
            data=node_df[node_df.original_type == show_label],
            size=label_size,
        )
        + theme_void()
    )
    return plot


def plot_kegg_global_map_plotnine(graph, hide_map=True):
    """Plot KEGG global map using plotnine.
    Parameters:
    -----------
    graph: KGML_graph class
        graph object of global map like ko01100, obtained by `KGML_graph()`.
    hide_map: bool
        if True, hide the map nodes.
    """
    node_df = graph.get_nodes()
    coords = graph.get_coords()
    if hide_map:
        node_df = node_df[node_df.original_type != "map"]
    plt = (
        ggplot()
        + geom_segment(
            aes(x="x", y="y", xend="xend", yend="yend"),
            color=coords["fgcolor"].tolist(),
            data=coords,
        )
        + geom_point(aes(x="x", y="y"), color=node_df["fgcolor"].tolist(), data=node_df)
        + theme_void()
    )
    return plt


def color_grad(
        minimum=-2, maximum=2, seq=0.1,
        min_col="#ffffff", max_col="#ff0000", round_num=2
):
    minmax = np.arange(minimum, maximum, seq)
    num_seq = len(minmax)
    conv = {}
    for num in range(num_seq):
        rounded = np.round(minmax[num], round_num)
        color1 = np.array(mpl.colors.to_rgb(min_col))
        color2 = np.array(mpl.colors.to_rgb(max_col))
        conv[rounded] = mpl.colors.to_hex(
            (1 - num / num_seq) * color1 + num / num_seq * color2
        )
    return conv


def hex2rgb(hex_str):
    """Convert hex string to rgb tuple.
    Parameters:
    ----------
    hex_str: str
        hex string, e.g. "#ffffff".
    """
    return tuple(int(hex_str.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
