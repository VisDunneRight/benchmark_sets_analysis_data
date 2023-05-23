import networkx as nx
import json
import os
import re


def read_storyline():
	for sname in os.listdir("../data/Storylines"):
		if os.path.splitext(sname)[1] == ".txt":
			print(sname)
			graph = {"nodes": [], "links": []}
			name_tsteps = {}
			with open(f"../data/Storylines/{sname}") as f:
				tstep = 0
				for ln in f.readlines():
					groups = re.split(r'[\t\n ]', ln.removesuffix('\n'))
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


if __name__ == '__main__':
	read_storyline()
