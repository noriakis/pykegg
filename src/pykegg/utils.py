import re
import os
import requests
import cv2
import warnings

import numpy as np
import matplotlib as mpl
import pandas as pd

from PIL import Image
from io import StringIO
from plotnine import (
    ggplot,
    geom_point,
    aes,
    geom_segment,
    theme_void,
    geom_rect,
    geom_text,
)

import pykegg


def overlay_opencv_image(
    node_df,
    path=None,
    pid=None,
    fill_color="color",
    transparent_colors=None,
    highlight_nodes=None,
    highlight_color="#ff0000",
    highlight_expand=2,
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
        the color number. Skip the node if `None`.
    transparent_color: list of str
        specify which color to be transparent.
        If `None`, default `["#FFFFFF", "#BFFFBF"]` is used.
    highlight_nodes: str
        the column in `node_df` specifying which nodes to be highlighted.
    hihglight_color: str
        the color of the highlighted nodes.
    highlight_expand: int
        the number of pixels to expand the highlighted nodes.
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
        if tmp_col is None:
            break
        if isinstance(tmp_col, list):
            num_col = len(tmp_col)
            nudge = tmp["width"] / num_col
            for col_num, one_tmp_col in enumerate(tmp_col):
                new_x0 = tmp["x0"] + (nudge * col_num)
                pt2 = (int(tmp["x"] + tmp["width"] / 2),
                       int(-1 * tmp["y"] + tmp["height"] / 2))
                canvas = cv2.rectangle(
                    img=canvas,
                    pt1=(
                        int(new_x0),
                        int(-1 * tmp["y0"])
                    ),
                    pt2 = pt2,
                    color=hex2rgb(one_tmp_col),
                    thickness=-1,
                )

        else:
            canvas = cv2.rectangle(
                img=canvas, rec=pos, color=hex2rgb(tmp_col), thickness=-1
            )

    ## Highlight the nodes by rectangle
    if highlight_nodes is not None:
        highlight_node_df = node_df[node_df[highlight_nodes]]
        for i in highlight_node_df.index:
            tmp = highlight_node_df.loc[i, :]
            pos = (
                int(tmp["x0"]),
                int(-1 * tmp["y0"]),
                int(tmp["width"]),
                int(tmp["height"]),
            )
            canvas = cv2.rectangle(
                img=canvas,
                pt1=(pos[0] - int(highlight_expand), pos[1] - int(highlight_expand)),
                pt2=(
                    pos[0] + pos[2] + int(highlight_expand),
                    pos[1] + pos[3] + int(highlight_expand),
                ),
                color=hex2rgb(highlight_color),
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
        node_df["graphics_name"] = node_df.graphics_name.apply(
            lambda x: x.split(",")[0]
        )

    ## Collapse subtypes
    edge_df_col = []
    for i in edge_df.index:
        tmp = edge_df.iloc[i, :]
        for subtype in tmp.subtypes:
            edge_df_col.append(
                [tmp.entry1, tmp.entry2, tmp.type, subtype, tmp.reaction]
            )
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
    low=-2, high=2, seq=0.01, low_col="#ffffff", high_col="#ff0000", round_num=2
):
    """Generate color gradient.

    Parameters:
    -----------
    low: float
        the lowest value.
    high: float
        the highest value.
    seq: float
        the sequence of the values.
    low_col: str
        the lowest color in HEX.
    high_col: str
        the highest color in HEX.
    round_num: int
        the number of digits to round.
    """
    minmax = np.arange(low, high, seq)
    num_seq = len(minmax)
    conv = {}
    for num in range(num_seq):
        rounded = np.round(minmax[num], round_num)
        color1 = np.array(mpl.colors.to_rgb(low_col))
        color2 = np.array(mpl.colors.to_rgb(high_col))
        conv[rounded] = mpl.colors.to_hex(
            (1 - (num + 1) / num_seq) * color1 + (num + 1) / num_seq * color2
        )
    return conv


def hex2rgb(hex_str):
    """Convert hex string to rgb tuple.

    Parameters:
    -----------
    hex_str: str
        hex string, e.g. "#ffffff".
    """
    return tuple(int(hex_str.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))


def deseq2_raw_map(
    results_df,
    path=None,
    pid=None,
    node_name_column="graphics_name",
    color_column="log2FoldChange",
    highlight_sig=False,
    highlight_color="#ff0000",
    highlight_padj_thresh=0.05,
):
    if ~highlight_sig:
        highlight_column = None

    sig_genes = results_df[results_df["padj"] < highlight_padj_thresh].index
    lfc_key = results_df[color_column].to_dict()

    graph = pykegg.KGML_graph(path=path, pid=pid)
    node_df = graph.get_nodes()

    if highlight_sig:
        highlight_value = []
        for node in node_df[node_name_column]:
            in_node = [i.replace("...", "") for i in node.split(",")]
            intersect = set(in_node) & set(sig_genes)
            if len(intersect) > 0:
                highlight_value.append(True)
            else:
                highlight_value.append(False)
        node_df["highlight"] = highlight_value
        highlight_column = "highlight"

    node_value = []

    for node in node_df[node_name_column]:
        ## Currently only graphics name is supported
        in_node = [i.replace("...", "") for i in node.split(",")]
        intersect = set(in_node) & set(lfc_key.keys())
        if len(intersect) > 0:
            tmp = [lfc_key[i] for i in lfc_key if i in intersect]
            node_value.append(np.mean(tmp))
        else:
            node_value.append(None)
    values = [n for n in node_value if n is not None]
    col_dic = color_grad2(
        low=min(values), mid=np.median(values), high=max(values), round_num=2, seq=0.01
    )

    node_df["color"] = [
        col_dic[np.round(x, 2)] if x is not None else None for x in node_value
    ]

    im_arr = Image.fromarray(
        overlay_opencv_image(
            node_df,
            path=path,
            pid=pid,
            highlight_nodes=highlight_column,
            highlight_color=highlight_color,
        )
    )
    return im_arr


def color_grad2(
    low=-2,
    mid=0,
    high=2,
    seq=0.01,
    low_col="#00ffff",
    mid_col="#ffffff",
    high_col="#ff0000",
    round_num=2,
):
    low_mid = np.arange(low, mid + seq, seq)
    mid_high = np.arange(mid, high + seq, seq)

    num_seq = len(low_mid)
    conv = {}
    for num in range(num_seq):
        rounded = np.round(low_mid[num], round_num)
        color1 = np.array(mpl.colors.to_rgb(low_col))
        color2 = np.array(mpl.colors.to_rgb(mid_col))
        conv[rounded] = mpl.colors.to_hex(
            (1 - (num + 1) / num_seq) * color1 + (num + 1) / num_seq * color2
        )
    num_seq = len(mid_high)
    for num in range(num_seq):
        rounded = np.round(mid_high[num], round_num)
        color1 = np.array(mpl.colors.to_rgb(mid_col))
        color2 = np.array(mpl.colors.to_rgb(high_col))
        conv[rounded] = mpl.colors.to_hex(
            (1 - (num + 1) / num_seq) * color1 + (num + 1) / num_seq * color2
        )
    return conv


def append_colors(
    node_df,
    candidate,
    new_column_name="color",
    candidate_column="graphics_name",
    delim=",",
    true_color="#ff0000",
    false_color="#ffffff",
):
    """Append colors to the node_df based on intersection with candidate ID list.

    Parameters:
    -----------
    node_df: DataFrame
        node data obtained by `get_nodes()`.
    candidate: list
        list of candidate IDs.
    new_column_name: str
        the name of the new column.
    candidate_column: str
        the column in `node_df` specifying candidate IDs.
    delim: str
        the delimiter of the node IDs. Typically "," for graphics_name,
        and " " for name.
    true_color: str
        the color of the candidate nodes.
    false_color: str
        the color of the non-candidate nodes.
    """
    if candidate_column=="name":
        delim=" "
    colors = []
    for i in node_df[candidate_column]:
        in_node = set([i.replace(" ", "").replace("...", "") for i in i.split(delim)])
        intersect = set(candidate) & in_node
        if len(intersect) > 0:
            colors.append(true_color)
        else:
            colors.append(false_color)
    node_df[new_column_name] = colors
    return node_df

def visualize_gseapy(gsea_res, colors,
                     pathway_name=None, pathway_id=None,
                     org="hsa",
                     column_name="graphics_name",
                     false_color="#707070"):
    if not isinstance(gsea_res, list):
        gsea_res = [gsea_res]
    if len(gsea_res) != len(colors):
        return
    
    ## Determine pathway name if not specified
    if pathway_name is None:
        if len(gsea_res) == 1:
            pathway_name = gsea_res[0].res2d.iloc[0,"Term"]
        else:
            pathway_names = []
            for res in gsea_res:
                pathway_names.append(set(res.res2d.Term))
            all_intersection = set.intersection(*pathway_names)
            if len(all_intersection) == 0:
                return
            first_res = gsea_res[0][gsea_res[0].Term.isin(list(all_intersection))]
            pathway_name = first_res.sort_values(by="P-value").Term.tolist()[0]
    
    ## [TODO] Fetch pathway ID given Term
    ## Still we need to specify organism name
    if pathway_id is None:
        pathway_id = pathway_name_to_id_dict(list_id=org)[pathway_name]
    
    graph = pykegg.KGML_graph(pid=pathway_id)
    nodes = graph.get_nodes()
    for node in nodes[column_name]:
        for e, res in enumerate(gsea_res):
            genes = res.res2d[res.res2d.Term==pathway_name].Genes.tolist()[0].split(";")
            nodes = append_colors(nodes, genes, new_column_name="color"+str(e),
                                  candidate_column=column_name,
                                  true_color=colors[e],
                                 false_color=false_color)
    
    col_col = [i for i in nodes.columns if i.startswith("color")]
    qc = list()
    for id in nodes.id:
        tmp = nodes[nodes.id==id].loc[:, col_col]
        tmp_colors = tmp.iloc[0,:].tolist()
        qc.append(tmp_colors)
    nodes["color"] = qc
    kegg_map = pykegg.overlay_opencv_image(nodes, pid=pathway_id)
    return kegg_map


def pathway_name_to_id_dict(list_id="hsa"):
    response = requests.get("https://rest.kegg.jp/list/pathway/"+list_id)
    check_cache(response)
    df = pd.read_csv(StringIO(
        response.content.decode("utf-8")
    ), sep="\t", header=None)
    df.index = df[1].apply(lambda x: x.split(" - ")[0])
    return df[0].to_dict()


def id_to_name_dict(list_id="hsa", column=3, semicolon=True, comma=True):
    response = requests.get("https://rest.kegg.jp/list/"+list_id)
    check_cache(response)

    df = pd.read_csv(StringIO(
        response.content.decode("utf-8")
    ), sep="\t", header=None)
    if semicolon:
        semicolon_df = df[3].apply(lambda x: x.split(";")[0])
    else:
        semicolon_df = df[3]
    if comma:
        comma_df = semicolon_df.apply(lambda x: x.split(",")[0])
    else:
        comma_df = semicolon_df
    comma_df.index = df[0]
    return comma_df.to_dict()


def check_cache(response):
    if not response.from_cache:
        warnings.warn("If it is not the first time fetching, please use requests_cache for caching")