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


def collection_distributions_4in1(path_to_folder):
	fpath = path_to_folder + "/clean"
	data = []
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			gr = read_json(fpath + "/" + gfile)
			dpt = {"Node Count": len(gr["nodes"]), "Edge Count": len(gr["links"])}
			degrees = {}
			for nd in gr["nodes"]:
				degrees[nd["id"]] = 0
			for edge in gr["links"]:
				for nd in edge["nodes"]:
					degrees[nd] += 1
			degrees = [v for v in degrees.values()]
			dpt["Mean Degree"] = sum(degrees) / len(degrees)
			dpt["Maximum Degree"] = int(max(degrees))
			data.append(dpt)
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
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


def construct_adjacency(json_graph):
	adj_list = {}
	for nd in json_graph["nodes"]:
		adj_list[nd["id"]] = []
	for ed in json_graph["links"]:
		for nd in ed["nodes"]:
			for nd_other in ed["nodes"]:
				if nd != nd_other:
					adj_list[nd].append(nd_other)
	return adj_list


def single_graph_charts(path_to_folder, binned=True):
	fpath = path_to_folder + "/clean"
	if "charts" not in os.listdir(path_to_folder):
		os.mkdir(path_to_folder + "/charts")
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			data = []
			gr = read_json(fpath + "/" + gfile)
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
			chart.save(path_to_folder + "/charts/two_in_one.svg", embed_options={'renderer': 'svg'})


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
				print("both dir/undir edges")
			if "clean_nx_json" not in os.listdir(path_to_folder):
				os.mkdir(path_to_folder + "/clean_nx_json")
			with open(path_to_folder + "/clean_nx_json/" + cfile, 'w') as fdesc:
				json.dump(new_g, fdesc, indent=2)


if __name__ == '__main__':
	collection = "../data/investment interdependence"
	# collection_distributions_4in1(collection)
	# single_graph_charts(collection, binned=False)
	# convert_all_nonmultigraph_to_nx_json(collection)

	# with open("../data/investment interdependence/clean_nx_json/investment_obstacles.json") as fds:
	# 	grz = json.load(fds)
	# gr = nx.node_link_graph(grz)
	# print(gr.nodes(data=True))
	# print(gr.edges(data=True))
	# print(gr.nodes[1]["obstacle_name"])
