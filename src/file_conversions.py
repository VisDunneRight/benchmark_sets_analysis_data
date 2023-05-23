import networkx as nx
import json
import os
import re


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
			with open(f"../data/Storylines/{sname.replace('.txt', '.json')}", 'w') as f:
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
			with open(f"../data/scotch/{sname.replace('.src', '.json')}", 'w') as f:
				json.dump(graph, f)


def read_rome():
	for gfile in os.listdir("../data/rome")[:1]:
		if os.path.splitext(gfile)[1] == ".graphml":
			g = nx.read_graphml(f"../data/rome/{gfile}")
			graph = {"nodes": [], "links": []}
			for v in g.nodes:
				graph["nodes"].append({"id": v.replace('n', '')})
			for e in g.edges:
				graph["links"].append({})


if __name__ == '__main__':
	# read_storyline()
	# read_scotch()
	read_rome()
