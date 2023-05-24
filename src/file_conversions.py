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


if __name__ == '__main__':
	# read_storyline()
	# read_scotch()
	# read_rome()
	# read_north()
	# read_kegg()
	# read_webcompute()
	read_airlines()


	# for fil in os.listdir("../data/north"):
	# 	if os.path.splitext(f"../data/north/{fil}")[1] == ".json":
	# 		shutil.move(f"../data/north/{fil}", f"../data/north/clean/{fil}")
