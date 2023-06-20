import os
import igraph
import requests

import pandas as pd
import numpy as np

from Bio.KEGG.KGML.KGML_parser import read

class KGML_graph:
    """KGML graph object."""
    def __init__(self, path=None, pid=None):
        """Initialize KGML graph object.

        Parameters:
        -----------
        path : str
            path to the KGML file.
        pid : str
            KEGG pathway identifier.
        """
        if path is not None and os.path.isfile(path):
            self.pathway = read(open(path, "r"))
        else:
            if pid is not None:
                # pid_pref = "".join(re.findall("[a-zA-Z]+", pid))
                path_cont = requests.get("https://rest.kegg.jp/get/" + pid + "/kgml")
                self.pathway = read(path_cont.content.decode())

    def get_graph(self, layout="native", add_group=False):
        """Get igraph object of the KGML graph.

        Parameters:
        -----------
        layout : str
            layout of the graph. If `native`, the original layout of the KGML file is used.
        """
        edge_df = self.get_edges(add_group=add_group)
        node_df = self.get_nodes()

        node_df["orig_id"] = node_df.id
        edge_df["entry1_orig_id"] = edge_df.entry1
        edge_df["entry2_orig_id"] = edge_df.entry2

        node_df.id = np.arange(0, node_df.shape[0], 1)
        node_df.index = node_df.id

        convert = dict(zip(node_df.orig_id, node_df.id))

        edge_df.entry1 = edge_df.entry1.replace(convert)
        edge_df.entry2 = edge_df.entry2.replace(convert)

        graph = igraph.Graph.DataFrame(edge_df, vertices=node_df)

        if layout != "native":
            layout = graph.layout(layout)
            layout = pd.DataFrame([l for l in layout])
            layout.columns = ["x", "y"]
            graph.vs["x"] = layout.x
            graph.vs["y"] = layout.y
        return graph

    def get_edges(self, add_group=False):
        """Get edges DataFrame of the KGML graph."""

        ## [TODO] Add edges indicating `in_group`
        rel_list = [
            [
                relation.entry1.id,
                relation.entry2.id,
                relation.type,
                relation.subtypes,
                None,
            ]
            for relation in self.pathway.relations
        ]
        if len(self.pathway.reactions) != 0:
            reacs = []
            for reaction in self.pathway.reactions:
                for substrate in reaction.substrates:
                    for product in reaction.products:
                        reacs.append(
                            [substrate.id, product.id, None, None, reaction.name]
                        )
            reacs = pd.DataFrame(reacs)
        else:
            reacs = None

        edges = pd.DataFrame(rel_list)
        edges = pd.concat([edges, reacs])

        if edges.shape[0] == 0:
            return None
        edges.columns = ["entry1", "entry2", "type", "subtypes", "reaction"]

        if add_group:
            group_df = []
            groups = [j for j in [self.pathway.entries[i] for i in self.pathway.entries] if j.type=="group"]
            if len(groups)==0:
                return edges
            for tmp_group in groups:
                for component in [i.id for i in tmp_group.components]:
                    group_df.append([tmp_group.id, component, "in_group", "in_group", None])
            group_df = pd.DataFrame(group_df)
            group_df.columns = ["entry1", "entry2", "type", "subtypes", "reaction"]

            edges = pd.concat([edges, group_df])
            
        edges.index = np.arange(0, edges.shape[0], 1)
        return edges

    def get_nodes(self, node_x_nudge=5, node_y_nudge=5,
        append_pathway_name=True):
        """Get nodes DataFrame of the KGML graph."""

        compounds = pd.DataFrame(
            [
                c[0]
                for c in [
                    [
                        [
                            j.id,
                            j.name,
                            element.x,
                            element.y,
                            j.name,
                            element.coords,
                            element.type,
                            element.width,
                            element.height,
                            element.fgcolor,
                            element.bgcolor,
                            [k.id for k in j.components],
                        ]
                        for element in j.graphics
                    ]
                    for j in [
                        self.pathway.entries[i]
                        for i in self.pathway.entries
                        if self.pathway.entries[i].type == "compound"
                    ]
                ]
            ]
        )
        compounds["original_type"] = "compound"

        if len(self.pathway.genes) != 0:
            genes = [
                [
                    [
                        gene.id,
                        gene.name,
                        element.x,
                        element.y,
                        element.name,
                        element.coords,
                        element.type,
                        element.width,
                        element.height,
                        element.fgcolor,
                        element.bgcolor,
                        None,
                    ]
                    for element in gene.graphics
                ]
                for gene in self.pathway.genes
            ]
            genes = pd.DataFrame([g[0] for g in genes])
            genes["original_type"] = "gene"
        else:
            genes = None

        if len(self.pathway.orthologs) != 0:
            orthos = [
                [
                    [
                        ortho.id,
                        ortho.name,
                        element.x,
                        element.y,
                        element.name,
                        element.coords,
                        element.type,
                        element.width,
                        element.height,
                        element.fgcolor,
                        element.bgcolor,
                        None,
                    ]
                    for element in ortho.graphics
                ]
                for ortho in self.pathway.orthologs
            ]
            orthos = pd.DataFrame([g[0] for g in orthos])
            orthos["original_type"] = "ortholog"
        else:
            orthos = None

        maps = [
            [
                [
                    m.id,
                    m.name,
                    element.x,
                    element.y,
                    element.name,
                    element.coords,
                    element.type,
                    element.width,
                    element.height,
                    element.fgcolor,
                    element.bgcolor,
                    None,
                ]
                for element in m.graphics
            ]
            for m in self.pathway.maps
        ]
        maps = pd.DataFrame([m[0] for m in maps])
        maps["original_type"] = "map"

        groups = [
            [
                [
                    j.id,
                    j.name,
                    element.x,
                    element.y,
                    j.name,
                    element.coords,
                    element.type,
                    element.width,
                    element.height,
                    element.fgcolor,
                    element.bgcolor,
                    [k.id for k in j.components],
                ]
                for element in j.graphics
            ]
            for j in [
                self.pathway.entries[i]
                for i in self.pathway.entries
                if self.pathway.entries[i].type == "group"
            ]
        ]
        groups = pd.DataFrame([g[0] for g in groups])
        groups["original_type"] = "group"

        nodes = pd.concat([genes, maps, groups, compounds, orthos])

        nodes.columns = [
            "id",
            "name",
            "x",
            "y",
            "graphics_name",
            "coords",
            "type",
            "width",
            "height",
            "fgcolor",
            "bgcolor",
            "group",
            "original_type",
        ]

        nodes["xmin"] = nodes["x"] - nodes["width"] + node_x_nudge
        nodes["xmax"] = nodes["x"] + nodes["width"] - node_x_nudge

        nodes["y"] = nodes["y"] * -1
        nodes["ymin"] = nodes["y"] - nodes["height"] + node_y_nudge
        nodes["ymax"] = nodes["y"] + nodes["height"] - node_y_nudge

        if append_pathway_name:
            nodes["pathway_name"] = self.pathway.name
        return nodes

    def get_coords(self):
        """Transform coords positions to edge DataFrame."""

        all_ortho = []
        for ortho in self.pathway.orthologs:
            graphics = ortho.graphics
            for element in graphics:
                seg_list = []
                if element.coords is None:
                    continue
                for i in range(len(element.coords)):
                    if i < len(element.coords) - 1:
                        seg_list.append(
                            (
                                element.coords[i][0],
                                element.coords[i][1],
                                element.coords[i + 1][0],
                                element.coords[i + 1][1],
                                element.name,
                                ortho.reaction,
                                element.type,
                                element.fgcolor,
                                element.bgcolor,
                            )
                        )
                all_ortho.append(pd.DataFrame(seg_list))
        if len(all_ortho) == 0:
            return None
        coords = pd.concat(all_ortho)
        coords.columns = [
            "x",
            "y",
            "xend",
            "yend",
            "name",
            "reaction",
            "type",
            "fgcolor",
            "bgcolor",
        ]
        coords["y"] = coords["y"] * -1
        coords["yend"] = coords["yend"] * -1
        return coords
