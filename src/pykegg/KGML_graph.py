import re
from Bio.KEGG.KGML.KGML_parser import read
import pandas as pd
import numpy as np
import igraph
from PIL import Image
import cv2
import matplotlib.pyplot as plt
import os
import requests

class KGML_graph:
    def __init__(self, path=None, pid=None):
        if path!=None and os.path.isfile(path):
            self.pathway = read(open(path, 'r'))
        else:
            if pid!=None:
                pid_pref = "".join(re.findall("[a-zA-Z]+", pid))
                path_cont = requests.get("https://rest.kegg.jp/get/"+pid+"/kgml")
                self.pathway = read(path_cont.content.decode())
                
    def get_graph(self, layout="native"):
        ed = self.get_edges()
        nd = self.get_nodes()
        
        nd["orig_id"] = nd.id
        ed["entry1_orig_id"] = ed.entry1
        ed["entry2_orig_id"] = ed.entry2
        
        nd.id = np.arange(0, nd.shape[0], 1)
        nd.index = nd.id
        
        convert = dict(zip(nd.orig_id, nd.id))
        
        ed.entry1 = ed.entry1.replace(convert)
        ed.entry2 = ed.entry2.replace(convert)
        
        g = igraph.Graph.DataFrame(ed, vertices=nd)
        
        if layout!="native":
            layout = g.layout(layout)
            layout = pd.DataFrame([l for l in layout])
            layout.columns = ["x", "y"]
            g.vs["x"] = layout.x
            g.vs["y"] = layout.y
        return(g)

    def get_edges(self):
        rel_list = [[r.entry1.id, r.entry2.id, r.type, r.subtypes, None] for r in self.pathway.relations]
        if (len(self.pathway.reactions)!=0):
            reacs = []
            for r in self.pathway.reactions:
                for s in r.substrates:
                    for p in r.products:
                        reacs.append([s.id, p.id, None, None, r.name])
            reacs = pd.DataFrame(reacs)
        else:
            reacs = None
            
        edges = pd.DataFrame(rel_list)
        edges = pd.concat([edges, reacs])
        edges.index = np.arange(0, edges.shape[0], 1)
        
        if (edges.shape[0]==0):
            return(None)
        edges.columns = ["entry1","entry2","type","subtypes","reaction"]
        
        return(edges)
    
    def get_nodes(self, node_x_nudge=5, node_y_nudge=5):

        
        compounds = pd.DataFrame([c[0] for c in [[[j.id, element.x, element.y, j.name, element.coords, element.type,
                        element.width, element.height, element.fgcolor, element.bgcolor,
            [k.id for k in j.components]] for element in j.graphics] 
           for j in [self.pathway.entries[i] for i in self.pathway.entries
                     if self.pathway.entries[i].type=="compound"]]])
        compounds["original_type"] = "compound"

        if len(self.pathway.genes)!=0:
            genes = [[[gene.id, element.x, element.y, element.name, element.coords, element.type,
                element.width, element.height, element.fgcolor, element.bgcolor, None] for element in gene.graphics] 
                     for gene in self.pathway.genes]
            genes = pd.DataFrame([g[0] for g in genes])
            genes["original_type"] = "gene"
        else:
            genes = None

        if len(self.pathway.orthologs)!=0:
            orthos = [[[ortho.id, element.x, element.y, element.name, element.coords, element.type,
                element.width, element.height, element.fgcolor, element.bgcolor, None] for element in ortho.graphics] 
                     for ortho in self.pathway.orthologs]
            orthos = pd.DataFrame([g[0] for g in orthos])
            orthos["original_type"] = "ortholog"
        else:
            orthos = None
            
        maps = [[[m.id, element.x, element.y, element.name, element.coords, element.type,
            element.width, element.height, element.fgcolor, element.bgcolor, None] for element in m.graphics] 
                 for m in self.pathway.maps]
        maps = pd.DataFrame([m[0] for m in maps])
        maps["original_type"] = "map"
        
        groups = [[[j.id, element.x, element.y, j.name, element.coords, element.type,
                                element.width, element.height, element.fgcolor, element.bgcolor,
                    [k.id for k in j.components]] for element in j.graphics] 
                   for j in [self.pathway.entries[i] for i in self.pathway.entries
                             if self.pathway.entries[i].type=="group"]]
        groups = pd.DataFrame([g[0] for g in groups])
        groups["original_type"] = "group"
        
        nodes = pd.concat([genes, maps, groups, compounds, orthos])
        
        nodes.columns = ["id","x","y","name","coords","type",
                        "width","height","fgcolor","bgcolor","group","original_type"]
        
        nodes["xmin"] = nodes["x"] - nodes["width"] + node_x_nudge
        nodes["xmax"] = nodes["x"] + nodes["width"] - node_x_nudge
        
        nodes["y"] = nodes["y"]*-1
        nodes["ymin"] = nodes["y"] - nodes["height"] + node_y_nudge
        nodes["ymax"] = nodes["y"] + nodes["height"] - node_y_nudge
        
        return(nodes)
    
    def get_coords(self):
        all_ortho=[]
        for ortho in self.pathway.orthologs:
            graphics = ortho.graphics
            for element in graphics:
                seg_list = []
                if element.coords==None:
                    break
                for i, co in enumerate(element.coords):
                    if i < len(element.coords)-1:
                        seg_list.append((
                            element.coords[i][0],
                            element.coords[i][1],
                            element.coords[i+1][0],
                            element.coords[i+1][1],
                            element.name,
                            element.type,
                            element.fgcolor,
                            element.bgcolor))
                all_ortho.append(pd.DataFrame(seg_list))
        if (len(all_ortho)==0):
            return(None)
        coords = pd.concat(all_ortho)
        coords.columns = ["x","y","xend","yend","name","type","fgcolor","bgcolor"]
        coords["y"] = coords["y"]*-1
        coords["yend"] = coords["yend"]*-1
        return(coords)