import json
import os
import altair as alt
import networkx as nx


def read_json(path_to_file):
	with open(path_to_file) as fp:
		gr_json = json.load(fp)
	return gr_json


def plot_histogram(data, x_axis, save_filepath):
	data_points = alt.Data(values=data)
	chart = alt.Chart(data_points).mark_bar().encode(
		x=alt.X(f"{x_axis}:Q").bin(),
		y=alt.Y('count()', axis=alt.Axis(title="Number of Occurrences in Collection"))
	)
	chart.save(save_filepath, embed_options={'renderer': 'svg'})


def collection_distributions(path_to_folder):
	fpath = path_to_folder + "/clean"
	node_counts = []
	edge_counts = []
	avg_degree_counts = []
	max_degree_counts = []
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			gr = read_json(fpath + "/" + gfile)
			node_counts.append({"Node Count": len(gr["nodes"])})
			edge_counts.append({"Edge Count": len(gr["links"])})
			degrees = {}
			for nd in gr["nodes"]:
				degrees[nd["id"]] = 0
			for edge in gr["links"]:
				for nd in edge["nodes"]:
					degrees[nd] += 1
			degrees = [v for v in degrees.values()]
			avg_degree_counts.append({"Mean Degree": sum(degrees) / len(degrees)})
			max_degree_counts.append({"Maximum Degree": int(max(degrees))})
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	plot_histogram(node_counts, "Node Count", path_to_folder + "/charts/node_counts.svg")
	plot_histogram(edge_counts, "Edge Count", path_to_folder + "/charts/edge_counts.svg")
	plot_histogram(avg_degree_counts, "Mean Degree", path_to_folder + "/charts/average_degree.svg")
	plot_histogram(max_degree_counts, "Maximum Degree", path_to_folder + "/charts/max_degree.svg")


def collection_distributions_4in1(path_to_folder, tree_coll=False):
	fpath = path_to_folder + "/clean_nx_json"
	data = []
	n_max_min = [1000000, 0]
	e_max_min = [1000000, 0]
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			gr = read_json(fpath + "/" + gfile)
			if len(gr["nodes"]) < n_max_min[0]:
				n_max_min[0] = len(gr["nodes"])
			if len(gr["nodes"]) > n_max_min[1]:
				n_max_min[1] = len(gr["nodes"])
			if len(gr["links"]) < e_max_min[0]:
				e_max_min[0] = len(gr["links"])
			if len(gr["links"]) > e_max_min[1]:
				e_max_min[1] = len(gr["links"])
			dpt = {"Node Count": len(gr["nodes"]), "Edge Count": len(gr["links"])}
			degrees = {}
			for nd in gr["nodes"]:
				degrees[nd["id"]] = 0
			for edge in gr["links"]:
				degrees[edge["source"]] += 1
				if "directed" not in gr or not gr["directed"]:
					degrees[edge["target"]] += 1
			degrees = [v for v in degrees.values()]
			dpt["Mean Degree"] = sum(degrees) / len(degrees)
			dpt["Maximum Degree"] = int(max(degrees))
			data.append(dpt)
	print(f"{n_max_min[0]}-{n_max_min[1]} nodes, {e_max_min[0]}-{e_max_min[1]} edges")
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	if tree_coll:
		base = alt.Chart().mark_bar().encode()
		chart = alt.hconcat(data=alt.Data(values=data))
		chart |= base.encode(x=alt.X(f"Node Count:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0F7B6C"))
		chart |= base.encode(x=alt.X(f"Edge Count:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0B6E99"))
		chart.save(path_to_folder + "/charts/four_in_one.svg", embed_options={'renderer': 'svg'})
	else:
		base = alt.Chart().mark_bar().encode()
		chart = alt.vconcat(data=alt.Data(values=data))
		r1 = alt.hconcat()
		r1 |= base.encode(x=alt.X(f"Node Count:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0F7B6C"))
		r1 |= base.encode(x=alt.X(f"Edge Count:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0B6E99"))
		chart &= r1
		r2 = alt.hconcat()
		r2 |= base.encode(x=alt.X(f"Mean Degree:Q").bin(), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#AD1A72"))
		r2 |= base.encode(x=alt.X(f"Maximum Degree:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#6940A5"))
		chart &= r2
		chart.save(path_to_folder + "/charts/four_in_one.svg", embed_options={'renderer': 'svg'})


def collection_distributions_4in1_swapped(path_to_folder, tree_coll=False):
	fpath = path_to_folder + "/clean_nx_json"
	data = []
	n_max_min = [1000000, 0]
	e_max_min = [1000000, 0]
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			gr = read_json(fpath + "/" + gfile)
			if len(gr["nodes"]) < n_max_min[0]:
				n_max_min[0] = len(gr["nodes"])
			if len(gr["nodes"]) > n_max_min[1]:
				n_max_min[1] = len(gr["nodes"])
			if len(gr["links"]) < e_max_min[0]:
				e_max_min[0] = len(gr["links"])
			if len(gr["links"]) > e_max_min[1]:
				e_max_min[1] = len(gr["links"])
			dpt = {"Node Count": len(gr["nodes"]), "Edge Count": len(gr["links"]), "File": gfile.removesuffix(".json")}
			degrees = {}
			for nd in gr["nodes"]:
				degrees[nd["id"]] = 0
			degrees[""] = 0
			for edge in gr["links"]:
				for nd in edge["nodes"]:
					degrees[nd] += 1
			degrees = [v for v in degrees.values()]
			dpt["Mean Degree"] = sum(degrees) / len(degrees)
			dpt["Maximum Degree"] = int(max(degrees))
			data.append(dpt)
	print(f"{n_max_min[0]}-{n_max_min[1]} nodes, {e_max_min[0]}-{e_max_min[1]} edges")
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	base = alt.Chart().mark_bar().encode().properties(width=300, height=300)
	chart = alt.vconcat(data=alt.Data(values=data))
	r1 = alt.hconcat()
	r1 |= base.encode(x=alt.X(f"File:O", axis=alt.Axis(title=None, labelAngle=-45)), y=alt.Y("Node Count:Q"), color=alt.value("#0F7B6C"))
	r1 |= base.encode(x=alt.X(f"File:O", axis=alt.Axis(title=None, labelAngle=-45)), y=alt.Y("Edge Count:Q"), color=alt.value("#0B6E99"))
	chart &= r1
	r2 = alt.hconcat()
	r2 |= base.encode(x=alt.X(f"File:O", axis=alt.Axis(title=None, labelAngle=-45)), y=alt.Y("Mean Degree:Q"), color=alt.value("#AD1A72"))
	r2 |= base.encode(x=alt.X(f"File:O", axis=alt.Axis(title=None, labelAngle=-45)), y=alt.Y("Maximum Degree:Q"), color=alt.value("#6940A5"))
	chart &= r2
	chart.save(path_to_folder + "/charts/four_in_one.svg", embed_options={'renderer': 'svg'})


def construct_adjacency(json_graph, multigraph=False):
	adj_list = {}
	for nd in json_graph["nodes"]:
		adj_list[nd["id"]] = []
	adj_list[''] = []
	if json_graph["directed"]:
		for ed in json_graph["links"]:
			adj_list[ed["source"]].append(ed["target"])
	else:
		for ed in json_graph["links"]:
			adj_list[ed["source"]].append(ed["target"])
			adj_list[ed["target"]].append(ed["source"])
	return adj_list


def single_graph_charts(path_to_folder, binned=True):
	fpath = path_to_folder + "/clean_nx_json"
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	ct = 0
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			data = []
			gr = read_json(fpath + "/" + gfile)
			print(f"{len(gr['nodes'])} nodes, {len(gr['links'])} edges")
			adj = construct_adjacency(gr)
			search_start = 0
			seen_nds = set()
			while len(seen_nds) != len(gr["nodes"]):
				compsize = 1
				i = search_start
				while gr["nodes"][i]["id"] in seen_nds:
					i += 1
				search_start = i
				dfsq = [gr["nodes"][i]["id"]]
				seen_nds.add(gr["nodes"][i]["id"])
				while dfsq:
					nd_next = dfsq.pop()
					for nd_adj in adj[nd_next]:
						if nd_adj not in seen_nds:
							seen_nds.add(nd_adj)
							dfsq.append(nd_adj)
							compsize += 1
				data.append({"Component Size": compsize, "Node Degree": -1})
			for nd in gr["nodes"]:
				data.append({"Component Size": -1, "Node Degree": len(adj[nd["id"]])})
			if binned:
				base = alt.Chart().mark_bar().encode()
				chart = alt.hconcat(data=alt.Data(values=data))
				chart |= base.encode(x=alt.X(f"Component Size:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0F7B6C")).transform_filter(alt.FieldEqualPredicate(field='Node Degree', equal=-1))
				chart |= base.encode(x=alt.X(f"Node Degree:Q").bin(minstep=1), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0B6E99")).transform_filter(alt.FieldEqualPredicate(field='Component Size', equal=-1))
			else:
				base = alt.Chart().mark_bar(size=20).encode()
				chart = alt.hconcat(data=alt.Data(values=data))
				chart |= base.encode(x=alt.X(f"Component Size:Q"), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0F7B6C")).transform_filter(alt.FieldEqualPredicate(field='Node Degree', equal=-1))
				chart |= base.encode(x=alt.X(f"Node Degree:Q"), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0B6E99")).transform_filter(alt.FieldEqualPredicate(field='Component Size', equal=-1))
			chart.save(path_to_folder + f"/charts/two_in_one{'' if ct == 0 else ct}.svg", embed_options={'renderer': 'svg'})
			ct += 1


def just_degree_distr(path_to_folder, binned=True, multigraph=False):
	fpath = path_to_folder + "/clean_nx_json"
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	ct = 0
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			data = []
			gr = read_json(fpath + "/" + gfile)
			print(f"{len(gr['nodes'])} nodes, {len(gr['links'])} edges")
			adj = construct_adjacency(gr, multigraph=multigraph)
			for nd in gr["nodes"]:
				data.append({"Node Degree": len(adj[nd["id"]])})
			if binned:
				chart = alt.Chart(data=alt.Data(values=data), title=gfile).mark_bar().encode(x=alt.X(f"Node Degree:Q").bin(minstep=1), y=alt.Y('count()'), color=alt.value("#0B6E99"))
			else:
				chart = alt.Chart(data=alt.Data(values=data)).mark_bar().encode(x=alt.X(f"Node Degree:Q"), y=alt.Y('count()'), color=alt.value("#0B6E99"))
			chart.save(path_to_folder + f"/charts/degree_distr{'' if ct == 0 else ct}.svg", embed_options={'renderer': 'svg'})
			ct += 1


def convert_all_nonmultigraph_to_nx_json(path_to_folder):
	for cfile in os.listdir(path_to_folder + "/clean"):
		if os.path.splitext(cfile)[1] == ".json":
			with open(path_to_folder + "/clean/" + cfile) as fdesc:
				gr = json.load(fdesc)
			new_g = {"multigraph": False, "directed": True, "nodes": gr["nodes"], "links": []}
			all_dir = True
			all_undir = True
			for link in gr["links"]:
				new_lk = {}
				for k, v in link.items():
					if k == "nodes":
						new_lk["source"] = v[0]
						new_lk["target"] = v[1]
					elif k == "directed":
						if v:
							all_undir = False
						else:
							all_dir = False
					else:
						new_lk[k] = v
				new_g["links"].append(new_lk)
			if all_dir:
				new_g["directed"] = True
			elif all_undir:
				new_g["directed"] = False
			else:
				print(f"{path_to_folder} both dir/undir edges")
			if "clean_nx_json" not in os.listdir(path_to_folder):
				os.mkdir(path_to_folder + "/clean_nx_json")
			with open(path_to_folder + "/clean_nx_json/" + cfile, 'w') as fdesc:
				json.dump(new_g, fdesc, indent=2)


def sparkline_data(path_to_folder, degree_distribution_instead=False):
	fpath = path_to_folder + "/clean_nx_json"
	data = []
	if not degree_distribution_instead:
		dataset_name = path_to_folder.replace("../data/", "")
		for gfile in os.listdir(fpath):
			if os.path.splitext(gfile)[1] == ".json":
				gr = read_json(fpath + "/" + gfile)
				data.append(len(gr["nodes"]))
	else:
		dataset_name = path_to_folder[:path_to_folder.index("clean") - 1]
		dataset_name = dataset_name.replace("../data/", "")
		gr = read_json(path_to_folder)
		dnodes = {}
		for ed in gr["links"]:
			if ed["source"] not in dnodes:
				dnodes[ed["source"]] = 0
			if ed["target"] not in dnodes:
				dnodes[ed["target"]] = 0
			dnodes[ed["source"]] += 1
			dnodes[ed["target"]] += 1
		for nd in gr["nodes"]:
			if nd["id"] not in dnodes:
				data.append(0)
		for v in dnodes.values():
			data.append(v)
	nbins = min(len(data), 25)
	nbins = 1 if nbins == 0 else nbins
	data.sort()
	maxsz = data[-1]
	minsz = data[0]
	stepsz = round(((maxsz - minsz + 13) / nbins) / 5) * 5 if (maxsz - minsz) / nbins > 3.5 else round((maxsz - minsz) / nbins)
	stepsz = round(((maxsz - minsz + 10) / nbins) / 50) * 50 if (maxsz - minsz) / nbins > 35 else stepsz
	stepsz = round(((maxsz - minsz + 100) / nbins) / 500) * 500 if (maxsz - minsz) / nbins > 350 else stepsz
	stepsz = round(((maxsz - minsz + 1000) / nbins) / 5000) * 5000 if (maxsz - minsz) / nbins > 3500 else stepsz
	stepsz = round(((maxsz - minsz + 10000) / nbins) / 50000) * 50000 if (maxsz - minsz) / nbins > 35000 else stepsz
	stepsz = 1 if stepsz == 0 else stepsz
	keybin = "num_nodes" if not degree_distribution_instead else "node_degree"
	data_out = {"min": minsz, "max": maxsz, "step_size": stepsz, "num_bins": 0, "bins": [], keybin: []}
	idx = 0
	i = 0
	while idx < len(data):
		n_g_in_bin = 0
		while idx < len(data) and data[idx] < (i + 1) * stepsz:
			n_g_in_bin += 1
			idx += 1
		data_out[keybin].append(n_g_in_bin)
		data_out["bins"].append(i * stepsz)
		data_out["num_bins"] += 1
		i += 1
	while data_out[keybin][-1] == 0:
		data_out[keybin].pop()
		data_out["bins"].pop()
		data_out["num_bins"] -= 1
	# while data_out[keybin][0] == 0:
	# 	data_out[keybin].pop(0)
	# 	data_out["bins"].pop(0)
	# 	data_out["num_bins"] -= 1
	print(dataset_name)
	print(data_out)


def convert_to_gexf(folder):
	foldpath = "../data/" + folder
	print(foldpath)
	ct = 0
	print("graphs created:", ct, end='')
	if os.path.isdir(foldpath + "/clean_nx_json"):
		for fl in os.listdir(foldpath + "/clean_nx_json"):
			with open(foldpath + '/clean_nx_json/' + fl, 'r') as fd1:
				grj = json.load(fd1)
			gr = nx.node_link_graph(grj)
			try:
				nx.write_gexf(gr, foldpath + "/" + folder + "_gexf/data/" + os.path.splitext(fl)[0] + ".gexf")
			except TypeError:
				pass
			ct += 1
			print('\b'*len(str(ct - 1)) + str(ct), end='')
	else:
		print("\bno graphs to create", end='')
	print('')


def convert_to_graphml(folder):
	foldpath = "../data/" + folder
	print(foldpath)
	ct = 0
	print("graphs created:", ct, end='')
	if os.path.isdir(foldpath + "/clean_nx_json"):
		for fl in os.listdir(foldpath + "/clean_nx_json"):
			with open(foldpath + '/clean_nx_json/' + fl, 'r') as fd1:
				grj = json.load(fd1)
			gr = nx.node_link_graph(grj)
			try:
				nx.write_graphml(gr, foldpath + "/" + folder + "_graphml/data/" + os.path.splitext(fl)[0] + ".graphml")
			except nx.exception.NetworkXError:
				pass
			ct += 1
			print('\b'*len(str(ct - 1)) + str(ct), end='')
	else:
		print("\bno graphs to create", end='')
	print('')


def convert_to_gml(folder):
	foldpath = "../data/" + folder
	print(foldpath)
	ct = 0
	print("graphs created:", ct, end='')
	if os.path.isdir(foldpath + "/clean_nx_json"):
		for fl in os.listdir(foldpath + "/clean_nx_json"):
			with open(foldpath + '/clean_nx_json/' + fl, 'r') as fd1:
				grj = json.load(fd1)
			gr = nx.node_link_graph(grj)
			try:
				nx.write_gml(gr, foldpath + "/" + folder + "_gml/data/" + os.path.splitext(fl)[0] + ".gml")
			except nx.exception.NetworkXError:
				pass
			ct += 1
			print('\b'*len(str(ct - 1)) + str(ct), end='')
	else:
		print("\bno graphs to create", end='')
	print('')


if __name__ == '__main__':
	collection = "../data/steinlib"
	# just_degree_distr(collection)
	sparkline_data(collection)
	# convert_all_nonmultigraph_to_nx_json(collection)
	# collection_distributions_4in1(collection)

	# for collection in os.listdir("../data"):
	# 	if os.path.isdir("../data/" + collection) and len(os.listdir("../data/" + collection + "/" + collection + "_gml/data")) == 0:
	# 		convert_to_gml(collection)

	# sp_datas = []
	# for fl in os.listdir("../data"):
	# 	if os.path.exists(f"../data/{fl}/clean"):
	# 		nfd = 0
	# 		ftt = ""
	# 		for fle in os.listdir(f"../data/{fl}/clean"):
	# 			if os.path.splitext(fle)[1] == ".json":
	# 				nfd += 1
	# 				ftt = f"../data/{fl}/clean/{fle}"
	# 		if nfd == 1:
	# 			sp_datas.append(sparkline_data(ftt, degree_distribution_instead=True))
	# 		else:
	# 			sp_datas.append(sparkline_data(f"../data/{fl}"))
	#
	# with open("sparkline_data.json", "w") as fd:
	# 	json.dump(sp_datas, fd, indent=2)

	# sparkline_data(collection, degree_distribution_instead=False)

	# collection_distributions_4in1_swapped(collection, tree_coll=False)
	# collection_distributions_4in1(collection, tree_coll=False)
	# single_graph_charts(collection, binned=False)
	# just_degree_distr(collection, binned=True, multigraph=False)
	# convert_all_nonmultigraph_to_nx_json(collection)

	# with open("../data/investment interdependence/clean_nx_json/investment_obstacles.json") as fds:
	# 	grz = json.load(fds)
	# gr = nx.node_link_graph(grz)
	# print(gr.nodes(data=True))
	# print(gr.edges(data=True))
	# print(gr.nodes[1]["obstacle_name"])
