import csv
import itertools
import networkx as nx
import json
import os
import re
import shutil
import xml.etree.ElementTree as ET
from Bio import Phylo
import newick
import graphviz
from scipy.io import loadmat
from scipy.sparse import find
import numpy as np


def read_storyline():
    for sname in os.listdir("../data/Storylines"):
        if os.path.splitext(sname)[1] == ".txt":
            graph = {"nodes": [], "links": []}
            name_tsteps = {}
            with open(f"../data/Storylines/{sname}") as f:
                tstep = 0
                for ln in f.readlines():
                    groups = re.split(r'[\t ]', ln.removesuffix('\n'))
                    group_names = [re.split(',', val) for val in groups]
                    for gidx, gp in enumerate(group_names):
                        for nm in gp:
                            graph["nodes"].append({"id": nm + "_" + str(tstep), "group": gidx + 1})
                            if nm in name_tsteps and name_tsteps[nm] == tstep - 1:
                                graph["links"].append(
                                    {"nodes": [nm + "_" + str(tstep - 1), nm + "_" + str(tstep)], "directed": True})
                            name_tsteps[nm] = tstep
                    tstep += 1
            with open(f"../data/Storylines/clean/{sname.replace('.txt', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_scotch():
    for sname in os.listdir("../data/scotch"):
        if os.path.splitext(sname)[1] == ".src":
            if os.path.splitext(sname)[0] + ".xyz" in os.listdir("../data/scotch"):
                positions = {}
                with open(f"../data/scotch/{sname.replace('src', 'xyz')}", 'r') as fx:
                    for ln in fx.readlines():
                        sep = re.split(r'[\t ]', ln.removesuffix('\n'))
                        if len(sep) > 1:
                            positions[sep[0]] = [v for v in sep[1:]]
                graph = {"nodes": [], "links": []}
                link_seen = set()
                with open(f"../data/scotch/{sname}") as f:
                    for ln in f.readlines():
                        sep = re.split(r'[\t ]', ln.removesuffix('\n'))
                        if len(sep) > 2:
                            ndict = {"id": sep[0], "value": sep[1], "x": positions[sep[0]][0],
                                     "y": positions[sep[0]][1]}
                            if len(positions[sep[0]]) == 3:
                                ndict["z"] = positions[sep[0]][2]
                            graph["nodes"].append(ndict)
                            for i in range(int(sep[2])):
                                lk = (sep[0], sep[2 * i + 4])
                                if (lk[1], lk[0]) not in link_seen:
                                    link_seen.add(lk)
                                    graph["links"].append(
                                        {"nodes": [lk[0], lk[1]], "value": sep[2 * i + 3], "directed": False})
            else:
                graph = {"nodes": [], "links": []}
                link_seen = set()
                with open(f"../data/scotch/{sname}") as f:
                    for ln in f.readlines():
                        sep = re.split(r'[\t ]', ln.removesuffix('\n'))
                        if len(sep) > 2:
                            graph["nodes"].append({"id": sep[0], "value": sep[1]})
                            for i in range(int(sep[2])):
                                lk = (sep[0], sep[2 * i + 4])
                                if (lk[1], lk[0]) not in link_seen:
                                    link_seen.add(lk)
                                    graph["links"].append(
                                        {"nodes": [lk[0], lk[1]], "value": sep[2 * i + 3], "directed": False})
            with open(f"../data/scotch/clean/{sname.replace('.src', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_rome():
    for gfile in os.listdir("../data/rome"):
        if os.path.splitext(gfile)[1] == ".graphml":
            g = nx.read_graphml(f"../data/rome/{gfile}")
            graph = {"nodes": [], "links": []}
            for v in g.nodes:
                graph["nodes"].append({"id": v.replace('n', '')})
            for e in g.edges:
                graph["links"].append({"nodes": [e[0].replace('n', ''), e[1].replace('n', '')], "directed": False})
            with open(f"../data/rome/clean/{gfile.replace('.graphml', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_north():
    for gfile in os.listdir("../data/north"):
        if os.path.splitext(gfile)[1] == ".graphml":
            g = nx.read_graphml(f"../data/north/{gfile}")
            graph = {"nodes": [], "links": []}
            for v in g.nodes:
                graph["nodes"].append({"id": v.replace('n', '')})
            for e in g.edges:
                graph["links"].append({"nodes": [e[0].replace('n', ''), e[1].replace('n', '')], "directed": True})
            with open(f"../data/north/clean/{gfile.replace('.graphml', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_kegg():
    for kfile in os.listdir("../data/KEGG pathways"):
        if os.path.splitext(kfile)[1] == ".xml":
            tree = ET.parse(f"../data/KEGG pathways/{kfile}")
            root = tree.getroot()
            graph = {"nodes": [], "links": []}
            for child in root:
                if child.tag == 'metabolites':
                    for metabolite in child:
                        name = metabolite.find('name').text
                        formula = metabolite.find('formula').text
                        description = metabolite.find('description').text
                        graph["nodes"].append({"id": name, "formula": formula, "desciption": description})
                elif child.tag == 'reactions':
                    for reaction in child:
                        id = reaction.find('id').text
                        name = reaction.find('name').text
                        reactant = reaction.find('reactant').text
                        product = reaction.find('product').text
                        reversible = False if reaction.find('reversible').text == "false" else True
                        subsystem = reaction.find('subsystem').text
                        graph["links"].append(
                            {"nodes": [reactant, product], "directed": True if reversible == "false" else False,
                             "reaction_id": id, "reaction_name": name, "reversible": reversible,
                             "subsystem": subsystem})
            with open(f"../data/KEGG pathways/clean/{kfile.replace('.xml', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_webcompute():
    for wfile in os.listdir("../data/webcompute"):
        if wfile != "clean":
            g = nx.read_gml("../data/webcompute/" + wfile, label=None)
            graph = {"nodes": [], "links": []}
            for v in g.nodes:
                graph["nodes"].append({"id": v})
            for e in g.edges:
                graph["links"].append({"nodes": [e[0], e[1]], "directed": True})
            with open(f"../data/webcompute/clean/{wfile.replace('.gml', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_airlines():
    for afile in os.listdir("../data/airlines-migration-air traffic"):
        if afile != "clean":
            with open(f"../data/airlines-migration-air traffic/{afile}") as f:
                gr = json.load(f)
                graph = {"nodes": [], "links": []}
                for nd in gr["nodes"]:
                    if "tooltip" in nd:
                        graph["nodes"].append({"id": nd["id"], "x": nd["x"], "y": nd["y"], "tooltip": nd["tooltip"]})
                    else:
                        graph["nodes"].append({"id": nd["id"], "x": nd["x"], "y": nd["y"]})
                for eg in gr["links"]:
                    graph["links"].append({"nodes": [eg["source"], eg["target"]], "directed": gr["directed"]})
            with open(f"../data/airlines-migration-air traffic/clean/{afile}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_chess():
    with open("../data/chess/chess.json") as f:
        gr = json.load(f)
        graph = {"nodes": [], "links": []}
        idct = 0
        for i, chessline in enumerate(gr):
            for j, move in enumerate(chessline):
                graph["nodes"].append({"id": idct, "piece": move["level"], "move_count": move["type"]})
                if j > 0:
                    graph["links"].append({"nodes": [idct - 1, idct], "directed": True})
                idct += 1
    with open(f"../data/chess/clean/chess.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_mid():
    with open("../data/MID/MIDB_5.0.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        graph = {"nodes": [], "links": []}
        cur_disp = 2
        last_disp_people = [[], []]
        last_disp_dates = ["", ""]
        seen_people = set()
        for ln in rdr:
            if ln[1] not in seen_people:
                graph["nodes"].append({"id": ln[1]})
                seen_people.add(ln[1])
            if int(ln[0]) != cur_disp:
                graph["links"].append(
                    {"nodes": [nd for subl in last_disp_people for nd in subl], "sideA": last_disp_people[0],
                     "sideB": last_disp_people[1], "directed": True, "start": last_disp_dates[0],
                     "end": last_disp_dates[1]})
                last_disp_people = [[], []]
            if int(ln[9]) == 1:
                last_disp_people[0].append(ln[1])
            else:
                last_disp_people[1].append(ln[1])
            cur_disp = int(ln[0])
            last_disp_dates = [f"{ln[4] if ln[4] != str(-9) else '?'}/{ln[3] if ln[3] != str(-9) else '?'}/{ln[5]}",
                               f"{ln[7] if ln[7] != str(-9) else '?'}/{ln[6] if ln[6] != str(-9) else '?'}/{ln[8]}"]
    with open(f"../data/MID/clean/midb.json", 'w') as f:
        json.dump(graph, f, indent=2)
    with open("../data/MID/MIDIP_5.01.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        graph = {"nodes": [], "links": []}
        cur_disp = 3551001
        last_disp_people = [[], []]
        last_disp_dates = ["", ""]
        seen_people = set()
        for ln in rdr:
            if ln[2] not in seen_people:
                graph["nodes"].append({"id": ln[2]})
                seen_people.add(ln[2])
            if int(ln[1]) != cur_disp:
                graph["links"].append(
                    {"nodes": [nd for subl in last_disp_people for nd in subl], "sideA": last_disp_people[0],
                     "sideB": last_disp_people[1], "directed": True, "start": last_disp_dates[0],
                     "end": last_disp_dates[1], "dispute_id": ln[0]})
                last_disp_people = [[], []]
            if int(ln[10]) == 1:
                last_disp_people[0].append(ln[2])
            else:
                last_disp_people[1].append(ln[2])
            cur_disp = int(ln[1])
            last_disp_dates = [f"{ln[5] if ln[5] != str(-9) else '?'}/{ln[4] if ln[4] != str(-9) else '?'}/{ln[6]}",
                               f"{ln[8] if ln[8] != str(-9) else '?'}/{ln[7] if ln[7] != str(-9) else '?'}/{ln[9]}"]
    with open(f"../data/MID/clean/midip.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_mid_storyline():
    with open("../data/MID/MIDB_5.0.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        graph = {"directed": True, "multigraph": False, "nodes": [], "links": []}
        cur_id = 1
        for ln in rdr:
            graph["nodes"].append({"id": cur_id, "country": ln[1], "disp_num": int(ln[0]), "time": f"{ln[4] if ln[4] != str(-9) else '?'}/{ln[3] if ln[3] != str(-9) else '?'}/{ln[5]}", "side_a": bool(int(ln[9]))})
            graph["nodes"].append({"id": cur_id + 1, "country": ln[1], "disp_num": int(ln[0]), "time": f"{ln[7] if ln[7] != str(-9) else '?'}/{ln[6] if ln[6] != str(-9) else '?'}/{ln[8]}", "side_a": bool(int(ln[9]))})
            graph["links"].append({"source": cur_id, "target": cur_id + 1, "disp_num": int(ln[0])})
            cur_id += 2
    with open(f"../data/MID/clean_nx_json/midb.json", 'w') as f:
        json.dump(graph, f, indent=2)
    with open("../data/MID/MIDIP_5.01.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        graph = {"directed": True, "multigraph": False, "nodes": [], "links": []}
        cur_id = 1
        country_last_seen = {}
        for ln in rdr:
            graph["nodes"].append({"id": cur_id, "country": ln[2], "disp_num": int(ln[0]), "ord_incidnum": int(ln[1]),
                                   "start": f"{ln[5] if ln[5] != str(-9) else '?'}/{ln[4] if ln[4] != str(-9) else '?'}/{ln[6]}",
                                   "end": f"{ln[8] if ln[8] != str(-9) else '?'}/{ln[7] if ln[7] != str(-9) else '?'}/{ln[9]}",
                                   "side_a": bool(int(ln[11]))})
            if ln[2] in country_last_seen and country_last_seen[ln[2]][1] == int(ln[0]):
                graph["links"].append({"source": country_last_seen[ln[2]][0], "target": cur_id, "disp_num": int(ln[0])})
            country_last_seen[ln[2]] = (cur_id, int(ln[0]))
            cur_id += 1
    with open(f"../data/MID/clean_nx_json/midip.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_greenhouse_gas():
    with open("../data/world-greenhouse-gas-emissions/ilp_case.json") as f:
        gr = json.load(f)
        graph = {"nodes": [], "links": []}
        for nd in gr["nodes"]:
            graph["nodes"].append({"id": nd["name"]})
        for ed in gr["links"]:
            graph["links"].append({"nodes": [ed["source"], ed["target"]], "directed": True, "value": ed["value"]})
    with open("../data/world-greenhouse-gas-emissions/clean/wri_data.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_tree_of_life():
    tree = ET.parse("../data/Tree of Life/mnVPRQ-xR-eU1dAsh6dWJA_phyloxml.xml")
    root = tree.getroot()
    graph = {"nodes": [], "links": []}
    found_names = set()

    def traverse(node):
        node_info = {'tag': node.tag, 'value': "", 'children': []}
        for child in node:
            if child.tag == "clade":
                node_info['children'].append(traverse(child))
            # cname = child.find("name")
            # if cname is not None:
            # 	graph["links"].append({"nodes": [node_info['value'], cname.text], "directed": True})
            elif child.tag == "name":
                xname = child.text
                while xname in found_names:
                    if not (xname[-8:] == "subclade" or "#" in xname):
                        xname += " subclade"
                    elif xname[-1] == 'e':
                        xname += '#2'
                    else:
                        nex = str(int(xname[xname.index('#') + 1:]) + 1)
                        xname = xname[:xname.index('#') + 1]
                        xname += nex
                found_names.add(xname)
                node_info["value"] = xname
                graph["nodes"].append({"id": xname})
        return node_info

    tree_structure = traverse(root)
    bfsq = tree_structure["children"]
    while bfsq:
        next_q = bfsq.copy()
        bfsq.clear()
        for nd in next_q:
            for cld in nd["children"]:
                graph["links"].append({"nodes": [nd['value'], cld["value"]], "directed": True})
                bfsq.append(cld)

    with open("../data/Tree of Life/clean/tree_of_life.json", 'w') as f:
        json.dump(graph, f, indent=2)


def create_complete_graphs():
    for i in range(5, 51):
        graph = {"nodes": [], "links": []}
        for j in range(i):
            graph["nodes"].append({"id": j})
        for j, k in itertools.combinations(list(range(i)), 2):
            graph["links"].append({"nodes": [j, k], "directed": False})
        with open(f"../data/complete graphs/clean/complete_K_{i}.json", 'w') as f:
            json.dump(graph, f, indent=2)


def create_complete_bipartite_graphs():
    for i in range(5, 41):
        for ip in range(5, 41):
            graph = {"nodes": [], "links": []}
            for j in range(i + ip):
                graph["nodes"].append({"id": j})
            for k in range(i):
                for kip in range(ip):
                    graph["links"].append({"nodes": [k, kip + i], "directed": False})
            with open(f"../data/complete bipartite graphs/clean/complete_K_{i}_{ip}.json", 'w') as f:
                json.dump(graph, f, indent=2)


def create_knowncr():
    # Ci x Cj
    for i in range(3, 8):
        j = i
        while i * j <= 250:
            graph = {"nodes": [], "links": []}
            for ip in range(i):
                for jp in range(j):
                    graph["nodes"].append({"id": f"{ip}_{jp}"})
            for ip in range(i):
                for jp in range(j):
                    graph["links"].append({"nodes": [f"{ip}_{jp}", f"{ip}_{(jp + 1) % j}"], "directed": False})
                    graph["links"].append({"nodes": [f"{ip}_{jp}", f"{(ip + 1) % i}_{jp}"], "directed": False})
            with open(f"../data/knownCR/clean/C{i}_x_C{j}.json", 'w') as f:
                json.dump(graph, f, indent=2)
            j += 1

    # Petersen graphs P(j, 2) and P(j, 3) for 9 <= j <= 125
    for i in range(2, 4):
        for j in range(9, 126):
            graph = {"nodes": [], "links": []}
            for k in range(j):
                graph["nodes"].append({"id": f"0_{k}"})
                graph["nodes"].append({"id": f"1_{k}"})
                graph["links"].append({"nodes": [f"0_{k}", f"0_{(k + 1) % j}"], "directed": False})
                graph["links"].append({"nodes": [f"0_{k}", f"1_{k}"], "directed": False})
                graph["links"].append({"nodes": [f"1_{k}", f"1_{(k + i) % j}"], "directed": False})
            with open(f"../data/knownCR/clean/P({j},{i}).json", 'w') as f:
                json.dump(graph, f, indent=2)

    # Gi x Pj
    for gfile in os.listdir("../data/knownCR"):
        if os.path.splitext(gfile)[1] == ".graph6":
            gr = nx.read_graph6(f"../data/knownCR/{gfile}")
            gr_num = gfile[gfile.index('_') + 1:gfile.index('.')]
            for j in range(3, 50):
                graph = {"nodes": [], "links": []}
                for ig in range(5):
                    for ij in range(j + 1):
                        graph["nodes"].append({"id": f"{ig}_{ij}"})
                for ig in range(5):
                    for ij in range(j):
                        graph["links"].append({"nodes": [f"{ig}_{ij}", f"{ig}_{ij + 1}"], "directed": False})
                for ij in range(j + 1):
                    for e1, e2 in gr.edges():
                        graph["links"].append({"nodes": [f"{e1}_{ij}", f"{e2}_{ij}"], "directed": False})
                with open(f"../data/knownCR/clean/G{gr_num}_x_P{j}.json", 'w') as f:
                    json.dump(graph, f, indent=2)

    # Gi x Cj
    for gfile in os.listdir("../data/knownCR"):
        if os.path.splitext(gfile)[1] == ".graph6":
            gr = nx.read_graph6(f"../data/knownCR/{gfile}")
            gr_num = gfile[gfile.index('_') + 1:gfile.index('.')]
            for j in range(3, 51):
                graph = {"nodes": [], "links": []}
                for ig in range(5):
                    for ij in range(j):
                        graph["nodes"].append({"id": f"{ig}_{ij}"})
                for ig in range(5):
                    for ij in range(j):
                        graph["links"].append({"nodes": [f"{ig}_{ij}", f"{ig}_{(ij + 1) % j}"], "directed": False})
                for ij in range(j):
                    for e1, e2 in gr.edges():
                        graph["links"].append({"nodes": [f"{e1}_{ij}", f"{e2}_{ij}"], "directed": False})
                with open(f"../data/knownCR/clean/G{gr_num}_x_C{j}.json", 'w') as f:
                    json.dump(graph, f, indent=2)


def read_evolution():
    for efile in os.listdir("../data/evolution"):
        if os.path.splitext(efile)[1] == ".nex":
            with open(f"../data/evolution/{efile}") as f:
                hostfound = False
                for ln in f.readlines():
                    if ln[:5] == "\tTREE":
                        if not hostfound:
                            hstr = "host"
                            hostfound = True
                        else:
                            hstr = "parasite"
                        treestr = ln[ln.index('=') + 2:-1]
                        graph = extract_tree_structure_newick(treestr)
                        with open(f"../data/evolution/clean/{os.path.splitext(efile)[0]}_{hstr}.json", 'w') as fd:
                            json.dump(graph, fd, indent=2)


def extract_tree_structure_newick(newick_string):
    tree = newick.loads(newick_string)
    graph = {"nodes": [], "links": []}

    def traverse(node):
        # Create a dictionary to store the node's information
        if node.name is None:
            node_info = {'name': f"root", 'children': []}
        else:
            node_info = {'name': node.name, 'children': []}
        for child in node.descendants:
            node_info['children'].append(traverse(child))
        return node_info

    tree_structure = traverse(tree[0])
    bfsq = [tree_structure]
    tree_structure["name"] = "root_0"
    d_ct = 1
    while bfsq:
        next_q = bfsq.copy()
        bfsq.clear()
        for nd in next_q:
            if "." in nd["name"]:
                graph["nodes"].append({"id": nd["name"], "host_node": nd["name"][:nd["name"].index(".")]})
            else:
                graph["nodes"].append({"id": nd["name"]})
            for cld in nd["children"]:
                if cld["name"] == 'root':
                    cld["name"] = f"root_{d_ct}"
                    d_ct += 1
                graph["links"].append({"nodes": [nd['name'], cld["name"]], "directed": True})
                bfsq.append(cld)
    return graph


def read_graphviz():
    for gfile in os.listdir("../data/graphviz examples"):
        if os.path.splitext(gfile)[1] == ".gv":
            print(gfile)
            gr = nx.node_link_data(nx.nx_pydot.read_dot("../data/graphviz examples/" + gfile))
            is_dir = gr["directed"]
            del gr["directed"]
            del gr["multigraph"]
            del gr["graph"]
            for ed in gr["links"]:
                ed["nodes"] = [ed["source"], ed["target"]]
                del ed["source"]
                del ed["target"]
                ed["directed"] = is_dir
            with open("../data/graphviz examples/clean/" + gfile.replace(".gv", ".json"), 'w') as fd:
                json.dump(gr, fd, indent=2)


def read_trade():
    with open("../data/Trade data/7bb3ede5-98b7-466d-a74a-f1a66d9de786_Data.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        gdp = {}
        name = {}
        for ln in rdr:
            if ln[1] != '':
                if ln[4] != '..':
                    gdp[ln[1]] = float(ln[4])
                name[ln[1]] = ln[0]
    with open("../data/Trade data/Trade_BEH0_1999_Export_2020Jan17.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        graph = {"nodes": [], "links": []}
        seen_nodes = set()
        trade_total = {}
        for line in rdr:
            if line[2] in name and line[9] in name:
                if line[2] not in seen_nodes:
                    if line[2] in gdp:
                        graph["nodes"].append({"id": line[2], "GDP_1999": gdp[line[2]], "country_name": name[line[2]]})
                    else:
                        graph["nodes"].append({"id": line[2], "country_name": name[line[2]]})
                    seen_nodes.add(line[2])
                if line[9] not in seen_nodes:
                    if line[9] in gdp:
                        graph["nodes"].append({"id": line[9], "GDP_1999": gdp[line[9]], "country_name": name[line[9]]})
                    else:
                        graph["nodes"].append({"id": line[9], "country_name": name[line[9]]})
                    seen_nodes.add(line[9])
                if (line[2], line[9]) not in trade_total and (line[9], line[2]) not in trade_total:
                    trade_total[(line[2], line[9])] = float(line[5])
                elif (line[9], line[2]) not in trade_total:
                    trade_total[(line[2], line[9])] += float(line[5])
                else:
                    trade_total[(line[9], line[2])] = float(line[5])
        for pair, val in trade_total.items():
            graph["links"].append({"nodes": [pair[0], pair[1]], "directed": False, "trade_value": val})
        with open("../data/Trade data/clean/trade_data.json", 'w') as fd:
            json.dump(graph, fd, indent=2)


def read_randdag():
    for gfile in os.listdir("../data/randDAG"):
        if os.path.splitext(gfile)[1] == ".graphml":
            g = nx.read_graphml(f"../data/randDAG/{gfile}")
            graph = {"nodes": [], "links": []}
            for v in g.nodes:
                graph["nodes"].append({"id": v.replace('n', '')})
            for e in g.edges:
                graph["links"].append({"nodes": [e[0].replace('n', ''), e[1].replace('n', '')], "directed": True})
            with open(f"../data/randDAG/clean/{gfile.replace('.graphml', '.json')}", 'w') as f:
                json.dump(graph, f, indent=2)


def read_investment():
    graph = {"nodes": [], "links": []}
    with open("../data/investment interdependence/malone_investment_obstacles.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        for ln in rdr:
            graph["nodes"].append({"id": int(ln[0]), "obstacle_name": ln[1]})
    with open("../data/investment interdependence/columbus_business_district_OSU.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        for ln in rdr:
            graph["links"].append({"nodes": [int(ln[0]), int(ln[1])], "directed": True})
    with open("../data/investment interdependence/clean/investment_obstacles.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_california():
    graph = {"nodes": [], "links": []}
    link_seen = set()
    with open("data\california\california.txt") as f:
        for ln in f.readlines():
            ln_lst = ln.split()
            if ln_lst[0] == 'n':
                graph["nodes"].append({"id": ln_lst[1], "url": ln_lst[2]})
            elif ln_lst[0] == 'e':
                link = (ln_lst[1], ln_lst[2])
                if link not in link_seen:
                    link_seen.add(link)
                    graph["links"].append({"nodes": [ln_lst[1], ln_lst[2]], "directed": False})
    with open(f"data\california\clean\california.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_collaborations():
    with open("../data/collaborations/IEEE VIS papers 1990-2021 - Main dataset.csv") as f:
        rdr = csv.reader(f)
        next(rdr)
        auth_seen = set()
        graph1 = {"nodes": [], "links": []}
        graph2 = {"nodes": [], "links": []}
        for ln in rdr:
            auths = re.split(r'[;|,]', ln[10])
            for i, auth in enumerate(auths):
                if auth not in auth_seen:
                    graph1["nodes"].append({"id": auth})
                    auth_seen.add(auth)
                for auth2 in auths[i + 1:]:
                    graph1["links"].append(
                        {"nodes": [auth, auth2], "directed": False, "paper_title": ln[2], "DOI": ln[3],
                         "conference": ln[0], "year": ln[1]})
            cites = re.split(r'[;|,]', ln[12])
            graph2["nodes"].append({"id": ln[3], "paper_title": ln[2], "conference": ln[0], "year": ln[1]})
            if ln[12] != "":
                for cite in cites:
                    graph2["links"].append({"nodes": [ln[3], cite], "directed": True})
        with open("../data/collaborations/clean/collaboration_network.json", 'w') as fd:
            json.dump(graph1, fd, indent=2)
        with open("../data/collaborations/clean/citation_network.json", 'w') as fd:
            json.dump(graph2, fd, indent=2)
    with open("../data/collaborations/Cpan_edge.csv") as f:
        rdr = csv.reader(f)
        seen = set()
        graph = {"nodes": [], "links": []}
        for ln in rdr:
            if int(ln[0]) not in seen:
                graph["nodes"].append({"id": int(ln[0])})
            if int(ln[1]) not in seen:
                graph["nodes"].append({"id": int(ln[1])})
            graph["links"].append({"nodes": [int(ln[0]), int(ln[1])], "directed": False})
        with open("../data/collaborations/clean/cpan_perl_module_users.json", 'w') as fd:
            json.dump(graph, fd, indent=2)


def read_codecommits():
    for cfile in os.listdir("../data/code/software"):
        with open("../data/code/software/" + cfile) as f:
            graph = {"nodes": [], "links": []}
            gr = json.load(f)
            for nd in gr["nodes"]:
                graph["nodes"].append(nd)
            for lk in gr["links"]:
                graph["links"].append({"nodes": [lk["source"], lk["target"]], "directed": True})
        with open("../data/code/clean/" + cfile, 'w') as fd:
            json.dump(graph, fd, indent=2)


def read_pi():
    seen_nodes = set()  # nodes are proteins

    graph = {"nodes": [], "links": []}
    with open("data\protein interactions\FriesCards.tsv") as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr)

        for ln in rdr:
            proteins = re.split(r';', ln[5])
            for protein in proteins:
                if protein not in seen_nodes:
                    seen_nodes.add(protein)
                    graph["nodes"].append({"id": protein})

            graph["links"].append({"nodes": [proteins[0], proteins[1]],
                                   "year": ln[1],
                                   "type": ln[2],
                                   "publication title": ln[3],
                                   "evidence": ln[4],
                                   "directed": True})

    with open("data\protein interactions\clean\protein_interactions_publications.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_blogs():
    files = ["corpus_ner_geo", "huffington", "wikinews"]
    categories = ["person", "location", "organization", "miscellaneous"]
    for file in files:
        graph = {"nodes": [], "links": []}
        seen_nodes = set()  # nodes are all topics across four categories

        with open(r"data\blogposts-tweets-forum\Blogposts\\" + file + ".tsv", encoding="utf-8") as f:
            rdr = csv.reader(f, delimiter="\t")
            next(rdr)

            for ln in rdr:
                topics_per_line = []

                for i in range(2, 5):  # columns corresponding to categories
                    topics_by_cat = re.split(r'\|', ln[i])
                    topics_per_line += topics_by_cat
                    for topic in topics_by_cat:
                        if topic not in seen_nodes:
                            seen_nodes.add(topic)
                            graph["nodes"].append({"id": topic,
                                                   "category": categories[i - 2]})

                # For space sake we save as a hypergraph as opposed to creating
                # cliques per line

                if file == files[2]:
                    source = "wikinews"
                else:
                    source = ln[0]

                graph["links"].append({
                    "nodes": topics_per_line,
                    "source": source,
                    "time": ln[1],
                    "directed": False,
                    "hyperedge": True
                })

                for i in range(len(topics_per_line)):
                    for j in range(i + 1, len(topics_per_line)):
                        if file == files[2]:
                            source = "wikinews"
                        else:
                            source = ln[0]
                        graph["links"].append({
                            "nodes": [topics_per_line[i], topics_per_line[j]],
                            "source": source,
                            "time": ln[1],
                            "directed": False
                        })

            with open(r"data\blogposts-tweets-forum\clean\Blogposts\{}.json".format(file), 'w') as f:
                json.dump(graph, f, indent=2)


def read_mooc():  # bipartite user - target
    graph = {"nodes": [], "links": []}
    nodes_seen = set()

    with open(r"data\blogposts-tweets-forum\MOOC\mooc_actions.tsv", encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr)

        for ln in rdr:
            user = "user" + ln[1]
            target = "target" + ln[2]

            for elem in [user, target]:
                if elem not in nodes_seen:
                    nodes_seen.add(elem)
                    graph["nodes"].append({"id": elem})

            graph["links"].append({
                "nodes": [user, target],
                "timestamp": ln[-1],
                "actionID": ln[0],
                "directed": False
            })

    with open(r"data\blogposts-tweets-forum\clean\MOOC\mooc.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_tweets():
    for gfile in os.listdir(r"data\blogposts-tweets-forum\tweets"):
        if os.path.splitext(gfile)[1] == ".wdnet":
            with open(r"data\blogposts-tweets-forum\tweets\\" + gfile) as f:
                graph = {"nodes": [], "links": []}
                nodes_seen = set()
                rdr = csv.reader(f, delimiter=" ")

                for ln in rdr:
                    for i, elem in enumerate(ln[1:-1]):  # elems consist of hashtags or user mentions
                        if elem not in nodes_seen:
                            nodes_seen.add(elem)
                            graph["nodes"].append({"id": elem})

                        for elem2 in ln[i + 2:-1]:
                            graph["links"].append({
                                "nodes": [elem, elem2],
                                "timestamp": ln[0],
                                "weight": ln[-1],
                                "directed": False
                            })

                with open(r"data\blogposts-tweets-forum\clean\tweets\\" + gfile.replace(".wdnet", ".json"), 'w') as f:
                    json.dump(graph, f, indent=2)

        with open(r"data\blogposts-tweets-forum\tweets\rugby tweets\pro12_mentions.csv", encoding="utf-8") as f:
            graph = {"nodes": [], "links": []}
            nodes_seen = set()

            rdr = csv.reader(f)
            next(rdr)

            for ln in rdr:
                for user in [ln[1], ln[2]]:
                    if user not in nodes_seen:
                        nodes_seen.add(user)
                        graph["nodes"].append({"id": user})

                graph["links"].append({"nodes": [ln[1], ln[2]],
                                       "timestamp": ln[0],
                                       "directed": True
                                       })

        with open(r"data\blogposts-tweets-forum\clean\tweets\\pro12_mentions.json", 'w') as f:
            json.dump(graph, f, indent=2)


def read_contacts():
    dir_path = "data\social network\Contacts"
    for gfile in os.listdir(dir_path):
        if os.path.splitext(gfile)[1] == ".txt":
            with open(dir_path + "\\" + gfile) as f:
                graph = {"nodes": [], "links": []}
                nodes_seen = set()
                rdr = csv.reader(f, delimiter="\t")
                next(rdr)

                for ln in rdr:
                    for user in [ln[1], ln[2]]:
                        if user not in nodes_seen:
                            nodes_seen.add(user)
                            graph["nodes"].append({"id": user})

                    graph["links"].append({"nodes": [ln[1], ln[2]],
                                           "timestamp": ln[0],
                                           "directed": False
                                           })

            with open("data\social network\clean\Contacts\\" + gfile.replace(".txt", ".json"), 'w') as f:
                json.dump(graph, f, indent=2)


def read_facebook100():
    dir_path = r"data\social network\facebook100\\"
    for gfile in os.listdir(dir_path):

        if os.path.splitext(gfile)[1] == ".mat" and not os.path.isfile(
                r"data\social network\clean\facebook100\\" + gfile.replace(".mat", ".json")):
            graph = {"nodes": [], "links": []}
            mat_data = loadmat(dir_path + gfile, squeeze_me=True)

            """
            "local_info" variable, one row per node: a student/faculty status
            flag, gender, major, second major/minor (if applicable), dorm/house,
            year, and high school. Missing data is coded 0.
            """
            all_features = mat_data["local_info"]
            adj = mat_data["A"]

            for i, n_att in enumerate(all_features):
                n_att = [int(n_att[x]) for x in range(len(n_att))]
                graph["nodes"].append({"id": i,
                                       "student/faculty flag": n_att[0],
                                       "gender": n_att[1],
                                       "major": n_att[2],
                                       "second major": n_att[3],
                                       "dorm": n_att[4],
                                       "year": n_att[5],
                                       "highschool ID": n_att[6]})

            row, col, val = find(adj)

            for i in range(len(row)):
                graph["links"].append({"nodes": [int(row[i]), int(col[i])],
                                       "directed": False
                                       })

            with open(r"data\social network\clean\facebook100\\" + gfile.replace(".mat", ".json"), 'w') as f:
                json.dump(graph, f, indent=2)


def read_vanDeBunt():
    dir_path = r"data\social network\vanDeBunt\\"
    nodes = []

    # VARS.DAT is a single file with node attr for all adjencencies over time
    with open(dir_path + "VARS.DAT") as f:
        """
        From Original data:
        gender (1 = F, 2 = M),
        program (2-year, 3-year, 4-year), 
        and smoking (1 = yes, 2 = no)
        """

        i = 0
        for ln in f.readlines():
            row = re.split("\s+", ln.strip())

            gender = {1: "F", 2: "M"}
            nodes.append({
                "id": i,
                "gender": gender[int(row[0])],
                "program (num years)": int(row[1]),
                "smoking": bool(-(int(row[2]) - 2))
            })

            i += 1
    """
    The networks are coded as 0 = unknown, 1 = best friend, 2 = friend, 3 = friendly relation, 4 = neutral,
    5 = troubled relation, 6 = item non-response, 9 = actor non-response. Note that 6 and 9 are missing data codes.
    """

    key_vals = {1: "bestfriend",
                2: "friend",
                3: "friendly relation",
                4: "neutral",
                5: "troubled relation",
                6: "item-non response",
                9: "actor non-response"}

    for gfile in os.listdir(dir_path):
        if os.path.splitext(gfile)[1] == ".DAT" and gfile != "VARS.DAT":
            graph = {"nodes": nodes, "links": []}
            with open(dir_path + gfile) as f:

                i = 0
                for ln in f.readlines():
                    # self reported data is directed, we remove unknowns from the network but keep 6, 9s as info may be valuable
                    # we do remove self loops whenever thry are 9
                    row = re.split("\s+", ln.strip())

                    for j in range(len(row)):
                        if not (row[i] == row[j] and int(row[j]) == 9):
                            if int(row[j]) != 0:
                                graph["links"].append({
                                    "nodes": [i, j],
                                    "realtionship": key_vals[int(row[j])],
                                    "directed": True
                                })
                    i += 1

            with open(r"data\social network\clean\vanDeBunt\\" + gfile.replace(".DAT", ".json"), 'w') as f:
                json.dump(graph, f, indent=2)


def add_karate():
    graph = {"nodes": [], "links": []}
    karate = nx.karate_club_graph()

    for nd in karate.nodes:
        graph["nodes"].append({"id": nd,
                               "club": karate.nodes[nd]["club"]})
    for u, v in karate.edges:
        # print(karate[u][v]) <- docs said there were weights, but I cannot seem to find them
        graph["links"].append({"nodes": [u, v],
                               "directed": False})

    with open(r"data\social network\clean\karate\karate.json", 'w') as f:
        json.dump(graph, f, indent=2)


def read_worldmaps():
    dir_path = r"data\world maps\\"
    subfolders = [r"NL\\", r"USA\\", r"WORLD\\"]

    for subfolder in subfolders:
        nodes = {}
        links = []
        for gfile in os.listdir(dir_path + subfolder):
            name, extension = os.path.splitext(gfile)
            file_id = re.split("_", name)[1]

            # get nodes and edges first
            if extension == ".tsv":
                if file_id in ["adjacencies", "locs"]:
                    with open(dir_path + subfolder + gfile) as f:
                        rdr = csv.reader(f, delimiter="\t")
                        header = next(rdr)

                        for ln in rdr:
                            if file_id == "adjacencies":
                                links.append({"nodes": [ln[0], ln[1]],
                                              "directed": False
                                              })
                            if file_id == "locs":
                                nodes[ln[0]] = {header[i]: ln[i] for i in range(1,len(ln))}

    # make graphs for each "Mixed" file also filling in other attributes from above
    for gfile in os.listdir(dir_path + subfolder):
        name, extension = os.path.splitext(gfile)
        file_id = re.split("_", name)[1]
        if extension == ".tsv" and file_id not in ["adjacencies", "locs"]:
            graph = {"nodes": [], "links": links}

            with open(dir_path + subfolder + gfile) as f:
                rdr = csv.reader(f, delimiter="\t")
                header = next(rdr)

                for ln in rdr:
                    nodes[ln[0]].update({header[i]: ln[i] for i in range(1,len(ln))})

                for k in nodes.keys(): #convert node dict to list of dicts per node
                    nodes[k].update({"id": k})
                graph["nodes"] += (list(nodes.values()))

            with open( r"data\world maps\clean\\"+ subfolder + gfile.replace(".tsv", ".json"), 'w') as f:
                json.dump(graph, f, indent=2)


def read_steinlib():
    for sfile in os.listdir("../data/steinlib"):
        if os.path.splitext(sfile)[1] == ".stp":
            with open("../data/steinlib/" + sfile) as fd:
                n_graph = {"directed": False, "multigraph": False, "graph": {}, "nodes": [], "links": []}
                status = 0
                for line in fd.readlines():
                    if "SECTION" in line or "Section" in line:
                        if "Graph" in line[8:]:
                            status = 1
                        elif "Terminals" in line[8:]:
                            status = 2
                        else:
                            status = 0
                    elif "END" in line or "End" in line:
                        status = 0
                    elif status == 1:
                        if "Edges" in line:
                            n_graph["directed"] = False
                        elif "Arcs" in line:
                            n_graph["directed"] = True
                        elif "Nodes" in line:
                            print(sfile)
                            n_nodes = int(line.split()[1])
                            for i in range(n_nodes):
                                n_graph["nodes"].append({"id": i+1, "terminal": False})
                        elif line[0] == "A":
                            arc = line.split()
                            n_graph["links"].append({"source": int(arc[1]), "target": int(arc[2]), "value": float(arc[3])})
                        elif line[0] == "E":
                            edge = line.split()
                            n_graph["links"].append({"source": int(edge[1]), "target": int(edge[2])})
                    elif status == 2:
                        if "Root" in line:
                            n_graph["graph"]["root"] = int(line.split()[1])
                        if line[0] == "T":
                            n_graph["nodes"][int(line.split()[1]) - 1]["terminal"] = True
            with open(f"../data/steinlib/clean_nx_json/{os.path.splitext(sfile)[0]}.json", 'w') as fd:
                json.dump(n_graph, fd, indent=2)
            with open("../data/steinlib/" + sfile) as fd:
                n_graph = {"nodes": [], "links": []}
                status = 0
                is_dir = True
                for line in fd.readlines():
                    if "SECTION" in line or "Section" in line:
                        if "Graph" in line[8:]:
                            status = 1
                        elif "Terminals" in line[8:]:
                            status = 2
                        else:
                            status = 0
                    elif "END" in line or "End" in line:
                        status = 0
                    elif status == 1:
                        if "Edges" in line:
                            is_dir = False
                        elif "Arcs" in line:
                            is_dir = True
                        elif "Nodes" in line:
                            print(sfile)
                            n_nodes = int(line.split()[1])
                            for i in range(n_nodes):
                                n_graph["nodes"].append({"id": i+1, "terminal": False})
                        elif line[0] == "A":
                            arc = line.split()
                            n_graph["links"].append({"nodes": [int(arc[1]), int(arc[2])], "directed": is_dir, "value": float(arc[3])})
                        elif line[0] == "E":
                            edge = line.split()
                            n_graph["links"].append({"nodes": [int(edge[1]), int(edge[2])], "directed": is_dir})
                    elif status == 2:
                        if line[0] == "T":
                            n_graph["nodes"][int(line.split()[1]) - 1]["terminal"] = True
            with open(f"../data/steinlib/clean/{os.path.splitext(sfile)[0]}.json", 'w') as fd:
                json.dump(n_graph, fd, indent=2)


def read_tol_big():
    n_graph = {"directed": True, "multigraph": False, "nodes": [], "links": []}
    with open("../data/Tree of Life/treeoflife_nodes.csv", 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        for ln in rdr:
            n_graph["nodes"].append({"id": int(ln[0]), "name": ln[1], "child_nodes": int(ln[2]), "leaf_node": int(ln[3]), "tolorg_link": int(ln[4]), "extinct": int(ln[5]), "confidence": int(ln[6]), "phylesis": int(ln[7])})
    with open("../data/Tree of Life/treeoflife_links.csv", 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        for ln in rdr:
            n_graph["links"].append({"source": int(ln[0]), "target": int(ln[1])})
    with open(f"../data/Tree of Life/clean_nx_json/tree_of_life.json", 'w') as fd:
        json.dump(n_graph, fd, indent=2)


def read_subways():
    for tfile in os.listdir("../data/subways/subway"):
        with open(f"../data/subways/subway/{tfile}", 'r') as fd:
            jsg = json.load(fd)
            with open(f"../data/subways/clean_nx_json/{tfile}", 'w') as fd2:
                json.dump(jsg, fd2, indent=2)


if __name__ == '__main__':
    # read_storyline()
    # read_scotch()
    # read_rome()
    # read_north()
    # read_kegg()
    # read_webcompute()
    # read_airlines()
    # read_chess()
    # read_mid()
    # read_greenhouse_gas()
    # read_tree_of_life()
    # create_complete_graphs()
    # create_complete_bipartite_graphs()
    # create_knowncr()
    # read_evolution()
    # read_graphviz()
    # read_trade()
    # read_investment()
    # read_randdag()
    # read_california()
    # read_collaborations()
    # read_codecommits()
    # read_pi()
    # read_blogs()
    # read_mooc()
    # read_tweets()
    # read_contacts()
    # read_facebook100()
    # read_vanDeBunt()
    # add_karate()
    # read_worldmaps()
    # read_steinlib()
    # read_tol_big()
    # read_mid_storyline()
    # read_subways()

    # direct = "../data/investment interdependence/clean"
    # for fle in os.listdir(direct):
    # 	if os.path.splitext(fle)[1] == ".json":
    # 		with open(direct + "/" + fle) as fd1:
    # 			gr = json.load(fd1)
    # 		with open(direct + "/" + fle, 'w') as fd1:
    # 			json.dump(gr, fd1, indent=2)

    exit()
