import csv
import itertools
import networkx as nx
import json
import os
import re
import shutil
import xml.etree.ElementTree as ET


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
								graph["links"].append({"nodes": [nm + "_" + str(tstep-1), nm + "_" + str(tstep)], "directed": True})
							name_tsteps[nm] = tstep
					tstep += 1
			with open(f"../data/Storylines/clean/{sname.replace('.txt', '.json')}", 'w') as f:
				json.dump(graph, f)


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
							ndict = {"id": sep[0], "value": sep[1], "x": positions[sep[0]][0], "y": positions[sep[0]][1]}
							if len(positions[sep[0]]) == 3:
								ndict["z"] = positions[sep[0]][2]
							graph["nodes"].append(ndict)
							for i in range(int(sep[2])):
								lk = (sep[0], sep[2*i + 4])
								if (lk[1], lk[0]) not in link_seen:
									link_seen.add(lk)
									graph["links"].append({"nodes": [lk[0], lk[1]], "value": sep[2*i + 3], "directed": False})
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
									graph["links"].append({"nodes": [lk[0], lk[1]], "value": sep[2 * i + 3], "directed": False})
			with open(f"../data/scotch/clean/{sname.replace('.src', '.json')}", 'w') as f:
				json.dump(graph, f)


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
				json.dump(graph, f)


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
				json.dump(graph, f)


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
						reversible = reaction.find('reversible').text
						subsystem = reaction.find('subsystem').text
						graph["links"].append({"nodes": [reactant, product], "directed": True if reversible == "false" else False, "reaction_id": id, "reaction_name": name, "reversible": reversible, "subsystem": subsystem})
			with open(f"../data/KEGG pathways/clean/{kfile.replace('.xml', '.json')}", 'w') as f:
				json.dump(graph, f)


def read_webcompute():
	for wfile in os.listdir("../data/webcompute"):
		if wfile != "clean":
			g = nx.read_gml("../data/webcompute/" + wfile, label=None)
			graph = {"nodes": [], "links": []}
			for v in g.nodes:
				graph["nodes"].append({"id": v})
			for e in g.edges:
				graph["links"].append({"nodes": [e[0], e[1]], "directed": True})
			with open(f"../data/webcompute/clean/{wfile.replace('.graphml', '.json')}", 'w') as f:
				json.dump(graph, f)


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
				json.dump(graph, f)


def read_chess():
	with open("../data/chess/chess.json") as f:
		gr = json.load(f)
		graph = {"nodes": [], "links": []}
		idct = 0
		for i, chessline in enumerate(gr):
			for j, move in enumerate(chessline):
				graph["nodes"].append({"id": idct, "piece": move["level"], "move_count": move["type"]})
				if j > 0:
					graph["links"].append({"nodes": [idct-1, idct], "directed": True})
				idct += 1
	with open(f"../data/chess/clean/chess.json", 'w') as f:
		json.dump(graph, f)


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
				graph["links"].append({"nodes": [nd for subl in last_disp_people for nd in subl], "sideA": last_disp_people[0], "sideB": last_disp_people[1], "directed": True, "start": last_disp_dates[0], "end": last_disp_dates[1]})
				last_disp_people = [[], []]
			if int(ln[9]) == 1:
				last_disp_people[0].append(ln[1])
			else:
				last_disp_people[1].append(ln[1])
			cur_disp = int(ln[0])
			last_disp_dates = [f"{ln[4] if ln[4] != str(-9) else '?'}/{ln[3] if ln[3] != str(-9) else '?'}/{ln[5]}", f"{ln[7] if ln[7] != str(-9) else '?'}/{ln[6] if ln[6] != str(-9) else '?'}/{ln[8]}"]
	with open(f"../data/MID/clean/midb.json", 'w') as f:
		json.dump(graph, f)
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
				graph["links"].append({"nodes": [nd for subl in last_disp_people for nd in subl], "sideA": last_disp_people[0], "sideB": last_disp_people[1], "directed": True, "start": last_disp_dates[0], "end": last_disp_dates[1], "dispute_id": ln[0]})
				last_disp_people = [[], []]
			if int(ln[10]) == 1:
				last_disp_people[0].append(ln[2])
			else:
				last_disp_people[1].append(ln[2])
			cur_disp = int(ln[1])
			last_disp_dates = [f"{ln[5] if ln[5] != str(-9) else '?'}/{ln[4] if ln[4] != str(-9) else '?'}/{ln[6]}", f"{ln[8] if ln[8] != str(-9) else '?'}/{ln[7] if ln[7] != str(-9) else '?'}/{ln[9]}"]
	with open(f"../data/MID/clean/midip.json", 'w') as f:
		json.dump(graph, f)


def read_greenhouse_gas():
	with open("../data/world-greenhouse-gas-emissions/ilp_case.json") as f:
		gr = json.load(f)
		graph = {"nodes": [], "links": []}
		for nd in gr["nodes"]:
			graph["nodes"].append({"id": nd["name"]})
		for ed in gr["links"]:
			graph["links"].append({"nodes": [ed["source"], ed["target"]], "directed": True, "value": ed["value"]})
	with open("../data/world-greenhouse-gas-emissions/clean/wri_data.json", 'w') as f:
		json.dump(graph, f)


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
						nex = str(int(xname[xname.index('#')+1:]) + 1)
						xname = xname[:xname.index('#')+1]
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
		json.dump(graph, f)


def create_complete_graphs():
	for i in range(5, 51):
		graph = {"nodes": [], "links": []}
		for j in range(i):
			graph["nodes"].append({"id": j})
		for j, k in itertools.combinations(list(range(i)), 2):
			graph["links"].append({"nodes": [j, k], "directed": False})
		with open(f"../data/complete graphs/clean/complete_{i}.json", 'w') as f:
			json.dump(graph, f)


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
	create_complete_graphs()


	# for fil in os.listdir("../data/north"):
	# 	if os.path.splitext(f"../data/north/{fil}")[1] == ".json":
	# 		shutil.move(f"../data/north/{fil}", f"../data/north/clean/{fil}")
