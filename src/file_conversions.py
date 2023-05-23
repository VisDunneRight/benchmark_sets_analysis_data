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
								graph["links"].append({"nodes": [nm + "_" + str(tstep-1), nm + "_" + str(tstep)]})
							name_tsteps[nm] = tstep
					tstep += 1
			with open(f"../data/Storylines/{sname.replace('.txt', '.json')}", 'w') as f:
				json.dump(graph, f)


def read_scotch():
	for sname in os.listdir("../data/scotch"):
		if os.path.splitext(sname)[1] == ".src":
			graph = {"nodes": [], "links": []}
			link_seen = set()
			with open(f"../data/scotch/{sname}") as f:
				for ln in f.readlines():
					sep = re.split(r'[\t ]', ln.removesuffix('\n'))
					if len(sep) > 2:
						graph["nodes"].append({"id": sep[0], "value": sep[1]})
						for i in range(int(sep[2])):
							lk = (sep[0], sep[2*i + 4])
							if (lk[1], lk[0]) not in link_seen:
								link_seen.add(lk)
								graph["links"].append({"nodes": [lk[0], lk[1]], "value": sep[2*i + 3]})
			with open(f"../data/scotch/{sname.replace('.src', '.json')}", 'w') as f:
				json.dump(graph, f)


if __name__ == '__main__':
	read_scotch()
