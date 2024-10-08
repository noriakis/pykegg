import re
import os
import warnings
import math
from io import StringIO

import requests
import cv2

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from PIL import Image

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
        If `None`, default `["#FFFFFF", "#BFFFBF", "#BFBFFF"]` is used.
    highlight_nodes: str
        the column in `node_df` specifying which nodes to be highlighted.
    hihglight_color: str
        the color of the highlighted nodes.
    highlight_expand: int
        the number of pixels to expand the highlighted nodes.
    """
    if transparent_colors is None:
        transparent_colors = ["#FFFFFF", "#BFFFBF", "#BFBFFF"]

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

    ## For making transparent image
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    canvas = np.zeros([image.shape[0], image.shape[1], 3], dtype=np.uint8)
    canvas.fill(255)

    dst = cv2.cvtColor(image, cv2.COLOR_BGR2BGRA)
    for col in transparent_colors:
        hex_str = col[1:7]
        cand_color = tuple(int(hex_str[i : i + 2], 16) for i in (0, 2, 4))
        mask = np.all(image[:, :, :] == list(cand_color), axis=-1)
        dst[mask, 3] = 0

    for i in node_df.id:
        tmp = node_df[node_df.id == i].iloc[0, :]
        
        if np.isnan(tmp.x0):
            continue
        pos = (
            int(tmp["x0"]),
            int(-1 * tmp["y0"]),
            int(tmp["width"]),
            int(tmp["height"]),
        )
        tmp_col = tmp[fill_color]
        if tmp_col is None:
            continue
        if isinstance(tmp_col, list):
            num_col = len(tmp_col)
            nudge = tmp["width"] / num_col
            for col_num, one_tmp_col in enumerate(tmp_col):
                new_x0 = tmp["x0"] + (nudge * col_num)
                pt2 = (
                    int(tmp["x"] + tmp["width"] / 2),
                    int(-1 * tmp["y"] + tmp["height"] / 2),
                )
                canvas = cv2.rectangle(
                    img=canvas,
                    pt1=(int(new_x0), int(-1 * tmp["y0"])),
                    pt2=pt2,
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
        for i in highlight_node_df.id:
            tmp = highlight_node_df[highlight_node_df.id == i]
            if np.isnan(tmp.x0):
                continue
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
        Image numpy array
    kegg_map: np.array
        Image numpy array
    """
    rects = cv2.cvtColor(rects, cv2.COLOR_BGR2RGB)
    rects = Image.fromarray(rects).convert("RGBA")
    kegg_map = cv2.cvtColor(kegg_map, cv2.COLOR_BGRA2RGBA)
    kegg_map = Image.fromarray(kegg_map).convert("RGBA")
    result_image = Image.alpha_composite(rects, kegg_map)
    return cv2.cvtColor(np.asarray(result_image), cv2.COLOR_RGBA2BGRA)


def return_segments(graph, node_df=None, edge_df=None):
    """Return edge dataframe to having xend and yend

    Parameters:
    -----------
    graph: KGML_graph
        KGML_graph class object
    node_df: DataFrame
        node data obtained by `get_nodes()`.
    edge_df: DataFrame
        edge data obtained by `get_edges()`.
    """
    if node_df is None:
        node_df = graph.get_nodes()
    if edge_df is None:
        edge_df = graph.get_edges()

    ## Subset to nodes
    in_nodes = node_df.id.tolist()
    edge_df = edge_df[
        edge_df.apply(lambda x: x.entry1 in in_nodes and x.entry2 in in_nodes, axis=1)
    ]

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
    return seg_df

def hex2rgb(hex_str):
    """Convert hex string to rgb tuple.

    Parameters:
    -----------
    hex_str: str
        hex string, e.g. "#ffffff".
    """
    return tuple(int(hex_str.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))


def return_color_bar(
    width=1,
    height=0.6,
    bottom=0.8,
    min_value=-2,
    max_value=2,
    two_slope=True,
    center_value=0,
    colors=None,
    label="Label",
):
    """Return color bar as a numpy array.

    Parameters:
    -----------
    width: int
        width of the color bar.
    height: int
        height of the color bar.
    bottom: int
        bottom of the color bar.
    min_value: int
        minimum value of the color bar.
    max_value: int
        maximum value of the color bar.
    two_slope: bool
        if True, use two slope norm.
    center_value: int
        center value of the color bar.
    colors: list
        list of colors.
    label: str
        label of the color bar.
    """
    mpl.use("Agg")
    fig, ax = plt.subplots(figsize=(width, height))
    fig.subplots_adjust(bottom=bottom)

    if colors is None:
        colors = ["#0000ff", "#ffffff", "#ff0000"]

    cmap_grad = mpl.colors.LinearSegmentedColormap.from_list("cmap_grad", colors)
    if two_slope:
        norm = mpl.colors.TwoSlopeNorm(
            vmin=min_value, vcenter=center_value, vmax=max_value
        )
    else:
        norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
    colbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap_grad),
        cax=ax,
        orientation="horizontal",
        label=label,
    )
    fig.canvas.draw()
    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    #     with_alpha = np.dstack((data,
    # 255*np.ones((legend.shape[0], legend.shape[1]), dtype=np.uint8)))
    mpl.pyplot.close()
    return data


def append_legend(
    image,
    min_value=-2,
    max_value=2,
    center_value=0,
    two_slope=True,
    colors=None,
    width=1,
    height=0.6,
    bottom=0.8,
    pos="topright",
    label="Label",
):
    """Add specified legend to image array

    Parameters:
    -----------
    image: numpy array
        Image numpy array
    min_value: float
        Minimum value of the color bar
    max_value: float
        Maximum value of the color bar
    center_value: float
        Center value of the color bar
    two_slope: bool
        If True, use two slope color bar
    colors: list
        List of colors, if None, use blue-white-red color bar
    width: float
        Width of the color bar
    height: float
        Height of the color bar
    bottom: float
        Bottom position of the color bar
    pos: str
        Position specification, "topright", "bottomright", "bottomleft", "topleft"
    label: str
        Label of the color bar
    """
    canvas = np.zeros([image.shape[0], image.shape[1], 3], dtype="uint8")
    canvas.fill(255)
    legend = return_color_bar(
        min_value=min_value,
        max_value=max_value,
        center_value=center_value,
        two_slope=two_slope,
        colors=colors,
        label=label,
        width=width,
        height=height,
        bottom=bottom,
    )

    if pos == "topright":
        ## topright
        canvas[
            0 : legend.shape[0],
            (canvas.shape[1] - legend.shape[1]) : canvas.shape[1],
            :,
        ] = legend
    elif pos == "bottomright":
        ## bottomright
        canvas[
            (canvas.shape[0] - legend.shape[0]) : canvas.shape[0],
            (canvas.shape[1] - legend.shape[1]) : canvas.shape[1],
            :,
        ] = legend
    elif pos == "bottomleft":
        ## bottomleft
        canvas[
            (canvas.shape[0] - legend.shape[0]) : canvas.shape[0],
            0 : legend.shape[1],
            :,
        ] = legend
    elif pos == "topleft":
        ## topleft
        canvas[0 : legend.shape[0], 0 : legend.shape[1], :] = legend
    else:
        return

    dst = cv2.cvtColor(canvas, cv2.COLOR_BGR2BGRA)
    for col in ["#FFFFFF"]:
        hex_str = col[1:7]
        cand_color = tuple(int(hex_str[i : i + 2], 16) for i in (0, 2, 4))
        mask = np.all(canvas[:, :, :] == list(cand_color), axis=-1)
        dst[mask, 3] = 0

    return overlay(image, dst)


def deseq2_raw_map(
    results_df,
    path=None,
    pid=None,
    node_name_column="graphics_name",
    delim=",",
    color_column="log2FoldChange",
    highlight_sig=False,
    highlight_color="#ff0000",
    highlight_padj_thresh=0.05,
    colors=None,
    two_slope=True,
    center_value=0,
    show_legend=True,
    legend_label=None,
    legend_position="topright",
    legend_width=1,
    legend_height=0.6,
    legend_bottom=0.8,
):
    """Plot PyDESeq2 results on KEGG pathway map

    Parameters:
    -----------
    results_df: pandas.DataFrame
        PyDESeq2 results dataframe
    path: str
        Path to the KEGG pathway map
    pid: str
        KEGG pathway ID
    node_name_column: str
        Column name of the node name
    delim: str
        Delimiter of the node name
    color_column: str
        Column name of the color
    highlight_sig: bool
        If True, highlight significant genes
    highlight_color: str
        Color of the highlight
    highlight_padj_thresh: float
        P-value threshold for the highlight
    colors: list
        List of colors, if None, use blue-white-red color bar
    two_slope: bool
        If True, use two slope color bar
    center_value: float
        Center value of the color bar
    show_legend: bool
        If True, show legend
    legend_label: str
        Label of the legend
    legend_position: str
        Position of the legend
    legend_width: float
        Width of the legend
    legend_height: float
        Height of the legend
    legend_bottom: float
        Bottom position of the legend
    """
    if ~highlight_sig:
        highlight_column = None

    sig_genes = results_df[results_df["padj"] < highlight_padj_thresh].index
    lfc_key = results_df[color_column].to_dict()

    graph = pykegg.KGML_graph(path=path, pid=pid)
    node_df = graph.get_nodes()

    if highlight_sig:
        highlight_value = []
        for node in node_df[node_name_column]:
            in_node = [i.replace("...", "") for i in node.split(delim)]
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
        in_node = [i.replace("...", "") for i in node.split(delim)]
        intersect = set(in_node) & set(lfc_key.keys())
        if len(intersect) > 0:
            tmp = [lfc_key[i] for i in lfc_key if i in intersect]
            node_value.append(np.mean(tmp))
        else:
            node_value.append(None)
    values = [n for n in node_value if n is not None]

    if colors is None:
        colors = ["#0000ff", "#ffffff", "#ff0000"]
    if center_value == "median":
        center_value = np.median(values)
    cmap_grad = mpl.colors.LinearSegmentedColormap.from_list("cmap_grad", colors)
    if two_slope:
        norm = mpl.colors.TwoSlopeNorm(
            vmin=min(values), vcenter=center_value, vmax=max(values)
        )
    else:
        norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
    node_df["color"] = [
        mpl.colors.to_hex(cmap_grad(norm(x))) if x is not None else None
        for x in node_value
    ]

    ## If mid-value are needed
    # col_dic = color_grad2(
    #     low=min(values), mid=np.median(values), high=max(values), round_num=2, seq=0.01
    # )

    # node_df["color"] = [
    #     col_dic[np.round(x, 2)] if x is not None else None for x in node_value
    # ]

    im_arr = overlay_opencv_image(
        node_df,
        path=path,
        pid=pid,
        highlight_nodes=highlight_column,
        highlight_color=highlight_color,
    )

    if legend_label is None:
        legend_label = color_column

    if show_legend:
        im_arr = append_legend(
            im_arr,
            min_value=min(values),
            max_value=max(values),
            colors=colors,
            center_value=center_value,
            pos=legend_position,
            two_slope=two_slope,
            width=legend_width,
            height=legend_height,
            bottom=legend_bottom,
            label=legend_label,
        )

    return Image.fromarray(im_arr)


def append_colors(
    node_df,
    candidate,
    new_column_name="color",
    candidate_column="graphics_name",
    delim=",",
    true_color="#ff0000",
    false_color="#ffffff",
):
    """Append discrete colors to the node_df based on intersection with candidate ID list.

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
    if candidate_column == "name":
        delim = " "
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


def append_colors_continuous_values(
    node_df,
    lfc_dict,
    node_name_column="graphics_name",
    new_color_column="color",
    delim=",",
    colors=None,
    two_slope=True,
    center_value="median",
    orig_value=None,
    fix_min=None,
    fix_max=None,
    fix_center=None
):
    """Append continuous colors to the node_df based on the values in dict.

    Parameters:
    -----------
    node_df: DataFrame
        node data obtained by `get_nodes()`.
    lfc_dict: dict
        dict of values.
    node_name_column: str
        the column in `node_df` specifying node IDs.
    new_color_column: str
        the name of the new column.
    delim: str
        the delimiter of the node IDs. Typically "," for graphics_name,
        and " " for name.
    colors: list
        the colors to be used. Default is ["#0000ff", "#ffffff", "#ff0000"].
    two_slope: bool
        whether to use two-slope color scheme. Default is True.
    center_value: str or float
        the center value of the color scheme. Default is "median".
    orig_value: str
        If specified, append the values used to compute color in the DataFrame.
    fix_min: float
        fixed minimum value to calculate color
    fix_max: float
        fixed maximum value to calculate color
    fix_center: float
        fixed center value to calculate color
        
    """
    node_value = []
    for node in node_df[node_name_column]:
        in_node = [i.replace("...", "") for i in node.split(delim)]
        intersect = set(in_node) & set(lfc_dict.keys())
        if len(intersect) > 0:
            tmp = [lfc_dict[i] for i in lfc_dict if i in intersect]
            node_value.append(np.nanmean(tmp))
        else:
            node_value.append(None)
    
    if orig_value is not None:
    	node_df[orig_value] = node_value
    
    values = [n for n in node_value if n is not None]

    if colors is None:
        colors = ["#0000ff", "#ffffff", "#ff0000"]

    if center_value == "median":
        center_value = np.nanmedian(values)

    if fix_min is not None:
        minval = fix_min
    else:
        minval = min(values)

    if fix_max is not None:
        maxval = fix_max
    else:
        maxval = max(values)

    if fix_center is not None:
        centval = fix_center
    else:
        centval = center_value

    cmap_grad = mpl.colors.LinearSegmentedColormap.from_list("cmap_grad", colors)
    if two_slope:
        norm = mpl.colors.TwoSlopeNorm(
            vmin=minval, vcenter=centval, vmax=maxval
        )
    else:
        norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    node_df[new_color_column] = [
        mpl.colors.to_hex(cmap_grad(norm(x))) if x is not None else None
        for x in node_value
    ]

    # col_dic = pykegg.color_grad2(
    #     low=min(values), mid=np.median(values), high=max(values), round_num=2, seq=0.01,
    #     low_col=low_col, mid_col=mid_col, high_col=high_col
    # )

    # node_df[new_color_column] = [
    #     col_dic[np.round(x, 2)] if x is not None else None for x in node_value
    # ]

    return node_df


def overlay_continuous_values_with_legend(
    node_df,
    value_dict,
    path=None,
    pid=None,
    node_name_column="graphics_name",
    delim=",",
    colors=None,
    legend_label="value",
    legend_position="topright",
    legend_width=1,
    legend_height=0.6,
    legend_bottom=0.8,
    transparent_colors=None,
    two_slope=True,
    center_value="median",
):
    """Obtain the raw image of pathway and color the nodes, return the overlaid image with legend.

    Parameters:
    -----------
    node_df: DataFrame
        node data obtained by `get_nodes()`.
    """
    if transparent_colors is None:
        transparent_colors = ["#FFFFFF", "#BFFFBF"]
    node_value = []
    for node in node_df[node_name_column]:
        ## Currently only graphics name is supported
        in_node = [i.replace("...", "") for i in node.split(delim)]
        intersect = set(in_node) & set(value_dict.keys())
        if len(intersect) > 0:
            tmp = [value_dict[i] for i in value_dict if i in intersect]
            node_value.append(np.mean(tmp))
        else:
            node_value.append(None)
    values = [n for n in node_value if n is not None]

    if colors is None:
        colors = ["#0000ff", "#ffffff", "#ff0000"]

    cmap_grad = mpl.colors.LinearSegmentedColormap.from_list("cmap_grad", colors)

    if center_value == "median":
        center_value = np.median(values)

    if two_slope:
        norm = mpl.colors.TwoSlopeNorm(
            vmin=min(values), vcenter=center_value, vmax=max(values)
        )
    else:
        norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))

    node_df["color"] = [
        mpl.colors.to_hex(cmap_grad(norm(x))) if x is not None else None
        for x in node_value
    ]

    if path is None and pid is None:
        ## Infer
        pid = (
            node_df["pathway_name"]
            .apply(lambda x: x.split(":")[1])
            .unique()
            .tolist()[0]
        )
    im_arr = overlay_opencv_image(
        node_df,
        path=path,
        pid=pid,
        transparent_colors=transparent_colors,
    )

    im_arr = append_legend(
        im_arr,
        min_value=min(values),
        max_value=max(values),
        colors=colors,
        pos=legend_position,
        center_value=center_value,
        width=legend_width,
        height=legend_height,
        bottom=legend_bottom,
        label=legend_label,
    )

    return Image.fromarray(im_arr)


def visualize(
    pathway_name,
    genes,
    db=None,
    org=None,
    column_name="graphics_name",
    false_color="#707070",
    true_color="#FA8072",
    output=None,
):
    """Output pathway image based on pathway *name* and gene symbol list

    Parameters:
    -----------
    pathway_name: str
        pathway name (not ID)
    genes: str
        list of genes
    db: str
        database name
    org: str
        if not specified db, the parameter will be used to convert pathway name to ID
    column_name: str
        column name to match for in node data
    true_color: str
        HEX specifying color for matched nodes
    false_color: str
        HEX specifying color for not matched nodes
    output: str
        output image file, default to None, meaning return the Image
    """
    if db is None and org is None:
        raise ValueError("Please specify db or org")

    if "Human" in db:
        org = "hsa"
    elif "Mouse" in db:
        org = "mmu"
    else:
        if org is None:
            return

    pathway_id = pykegg.pathway_name_to_id_dict(list_id=org)[pathway_name]

    graph = pykegg.KGML_graph(pid=pathway_id)
    nodes = graph.get_nodes()
    nodes = nodes[nodes.original_type == "gene"]

    ## Append color to node dataframe
    nodes = append_colors(
        nodes,
        genes,
        true_color=true_color,
        false_color=false_color,
    )

    kegg_map_image = Image.fromarray(pykegg.overlay_opencv_image(nodes, pid=pathway_id))
    if output is not None:
        kegg_map_image.save(output)
        return
    else:
        return kegg_map_image


def visualize_gseapy(
    gsea_res,
    colors,
    pathway_name=None,
    pathway_id=None,
    org="hsa",
    column_name="graphics_name",
    false_color="#707070",
):
    """Visualize GSEApy results.

    Parameters:
    -----------
    gsea_res: GSEApy object or list of GSEApy objects
        GSEApy results.
    colors: str or list of str
        Colors to use for each gsea results.
    pathway_name: str
        Pathway name.
    pathway_id: str
        Pathway ID.
    org: str
        KEGG organism ID.
    column_name: str
        Column name to use for visualization.
    false_color: str
        Color to use for false nodes.
    """
    if not isinstance(gsea_res, list):
        gsea_res = [gsea_res]
    if not isinstance(colors, list):
        colors = [colors]
    if len(gsea_res) != len(colors):
        return

    ## Determine pathway name if not specified
    if pathway_name is None:
        if len(gsea_res) == 1:
            pathway_name = gsea_res[0].res2d.iloc[0, "Term"]
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
        pathway_id = pykegg.pathway_name_to_id_dict(list_id=org)[pathway_name]

    graph = pykegg.KGML_graph(pid=pathway_id)
    nodes = graph.get_nodes()
    for node in nodes[column_name]:
        for e, res in enumerate(gsea_res):
            genes = (
                res.res2d[res.res2d.Term == pathway_name].Genes.tolist()[0].split(";")
            )
            nodes = append_colors(
                nodes,
                genes,
                new_column_name="color" + str(e),
                candidate_column=column_name,
                true_color=colors[e],
                false_color=false_color,
            )

    col_col = [i for i in nodes.columns if i.startswith("color")]
    qc = []
    for node_id in nodes.id:
        tmp = nodes[nodes.id == node_id].loc[:, col_col]
        tmp_colors = tmp.iloc[0, :].tolist()
        qc.append(tmp_colors)
    nodes["color"] = qc
    kegg_map = pykegg.overlay_opencv_image(nodes, pid=pathway_id)
    return kegg_map


def pathway_name_to_id_dict(list_id="hsa"):
    """Get pathway name to ID dictionary.

    Parameters:
    -----------
    list_id: str
        organism ID.
    """
    response = requests.get("https://rest.kegg.jp/list/pathway/" + list_id)
    check_cache(response)
    df = pd.read_csv(StringIO(response.content.decode("utf-8")), sep="\t", header=None)
    df.index = df[1].apply(lambda x: x.split(" - ")[0])
    return df[0].to_dict()


def id_to_name_dict(list_id="hsa", column=3, semicolon=True, comma=True):
    """Get KEGG ID to name dictionary.

    Parameters:
    -----------
    list_id: str
        organism ID.
    column: int
        Column to use for name.
    semicolon: bool
        Whether to split by semicolon.
    comma: bool
        Whether to split by comma.
    """
    response = requests.get("https://rest.kegg.jp/list/" + list_id)
    check_cache(response)

    df = pd.read_csv(StringIO(response.content.decode("utf-8")), sep="\t", header=None)
    if semicolon:
        semicolon_df = df[column].apply(lambda x: x.split(";")[0])
    else:
        semicolon_df = df[column]
    if comma:
        comma_df = semicolon_df.apply(lambda x: x.split(",")[0])
    else:
        comma_df = semicolon_df
    comma_df.index = df[0]
    return comma_df.to_dict()


def check_cache(response):
    """Check if response is from cache.

    Parameters:
    -----------
    response: requests.Response
        Response object.
    """
    try:
        if not response.from_cache:
            warnings.warn(
                "If it is not the first time fetching, please use requests_cache for caching"
            )
    except:
        warnings.warn(
            "If it is not the first time fetching, please use requests_cache for caching and run install_cache()"
        )


def shorten_end(row, pct=0.8, absolute=None):
    """shorten segments by moving yend and xend


    Parameters:
    -----------
    row: pd.Series
        Series of edge data frame
    pct: float
        scaling factor
    absolute: float
        absolute distance to shorten
    """
    radians = math.atan2(row["yend"] - row["y"], row["xend"] - row["x"])
    dists = math.dist([row["x"], row["y"]], [row["xend"], row["yend"]])
    if absolute is None:
        short = dists * pct
    else:
        short = dists - absolute

    new_xend = row["x"] + short * math.cos(radians)
    new_yend = row["y"] + short * math.sin(radians)

    row["xend"] = new_xend
    row["yend"] = new_yend

    return row


def convert_id(x, c_dic, first_only=True):
    """convert ID based on dict

    Parameters:
    -----------
    x: str
        node name
    c_dic: dict
        dictionary obtained typically obtained by `id_to_name_dict`
        keys correspond to KEGG ID and values correspond to name
    first_only:
        return only first string separated by space

    """

    in_node = x.split(" ")
    if first_only:
        if in_node[0] in c_dic.keys():
            converted = c_dic[in_node[0]]
        else:
            converted = np.nan
    else:
        tmp_node = []
        for node in in_node:
            if node in c_dic.keys():
                tmp_node.append(c_dic[node])
            else:
                pass
        if len(tmp_node) != 0:
            converted = " ".join(tmp_node)
        else:
            converted = np.nan
    return converted


def parallel_edges(df, move_param=5):
    """Experimental function moving x and y positions if multiple edges are to be plotted in plotnine based on whether the y position is the same between two points

    Parameters
    ----------
    df: pd.DataFrame
        data frame returned by `return_segments`
    move_params: float or int
        parameter to control edge nudge
    """
    ## Identify multiple edges
    df_dup = df[df.subtypes.apply(lambda x: len(x)) > 1]
    df_dup_collapse = []
    new_col = df_dup.columns
    for row in df_dup.index:
        tmp_group = df_dup.loc[row, :]
        for i in tmp_group["subtypes"]:
            tmp_group_rep = tmp_group.copy()
            tmp_group_rep["subtypes"] = [i]
            df_dup_collapse.append(tmp_group_rep)
    df_dup_collapse = pd.concat(df_dup_collapse, axis=1).T
    df_nodup = df[df.subtypes.apply(lambda x: len(x)) == 1]

    df_dup_re = []
    for group in df_dup_collapse.groupby(["entry1", "entry2"]):
        tmp_group = group[1]
        radians = math.atan2(
            tmp_group.yend.unique() - tmp_group.y.unique(),
            tmp_group.xend.unique() - tmp_group.x.unique(),
        )
        nudge = np.linspace(-1 * move_param, move_param, tmp_group.shape[0])
        if tmp_group.y.unique()[0] == tmp_group.yend.unique()[0]:
            tmp_group["y"] = tmp_group["y"] + nudge
            tmp_group["yend"] = tmp_group["yend"] + nudge
        else:
            tmp_group["x"] = tmp_group["x"] + nudge
            tmp_group["xend"] = tmp_group["xend"] + nudge
        df_dup_re.append(tmp_group)
    seg_df = pd.concat([df_nodup, pd.concat(df_dup_re)])
    for change in ["xend", "x", "yend", "y"]:
        seg_df[change] = seg_df[change].apply(lambda x: float(x))

    return seg_df


def parallel_edges2(df, move_param=5):
    """Experimental function moving x and y positions if multiple edges are to be plotted in plotnine based on degrees between points

    Parameters
    ----------
    df: pd.DataFrame
        data frame returned by `return_segments`
    move_params: float or int
        parameter to control edge nudge
    """
    ## Identify multiple edges
    df_dup = df[df.subtypes.apply(lambda x: len(x)) > 1]
    df_dup_collapse = []
    new_col = df_dup.columns
    for row in df_dup.index:
        tmp_group = df_dup.loc[row, :]
        for i in tmp_group["subtypes"]:
            tmp_group_rep = tmp_group.copy()
            tmp_group_rep["subtypes"] = [i]
            df_dup_collapse.append(tmp_group_rep)
    df_dup_collapse = pd.concat(df_dup_collapse, axis=1).T
    df_nodup = df[df.subtypes.apply(lambda x: len(x)) == 1]

    df_dup_re = []
    for group in df_dup_collapse.groupby(["entry1", "entry2"]):
        tmp_group = group[1]
        abs_deg = abs(
            math.degrees(
                math.atan2(
                    tmp_group.yend.unique() - tmp_group.y.unique(),
                    tmp_group.xend.unique() - tmp_group.x.unique(),
                )
            )
        )
        nudge = np.linspace(-1 * move_param, move_param, tmp_group.shape[0])
        if abs_deg < 45:
            tmp_group["y"] = tmp_group["y"] + nudge
            tmp_group["yend"] = tmp_group["yend"] + nudge
        else:
            tmp_group["x"] = tmp_group["x"] + nudge
            tmp_group["xend"] = tmp_group["xend"] + nudge
        df_dup_re.append(tmp_group)
    seg_df = pd.concat([df_nodup, pd.concat(df_dup_re)])
    for change in ["xend", "x", "yend", "y"]:
        seg_df[change] = seg_df[change].apply(lambda x: float(x))

    return seg_df
