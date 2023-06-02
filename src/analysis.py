import json
import os
import altair as alt


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
	r1 |= base.encode(x=alt.X(f"Node Count:Q").bin(), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0F7B6C"))
	r1 |= base.encode(x=alt.X(f"Edge Count:Q").bin(), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#0B6E99"))
	chart &= r1
	r2 = alt.hconcat()
	r2 |= base.encode(x=alt.X(f"Mean Degree:Q").bin(), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#AD1A72"))
	r2 |= base.encode(x=alt.X(f"Maximum Degree:Q").bin(), y=alt.Y('count()', axis=alt.Axis(title=None)), color=alt.value("#6940A5"))
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


def single_graph_charts(path_to_folder):
	fpath = path_to_folder + "/clean"
	data = []
	for gfile in os.listdir(fpath):
		if os.path.splitext(gfile)[1] == ".json":
			comp_sizes = []
			gr = read_json(fpath + "/" + gfile)
			adj = construct_adjacency(gr)
			seen_nds = set()
			while len(seen_nds) != len(gr["nodes"]):
				dfsq = [gr["nodes"][0]["id"]]
				while dfsq:
					nd_next = 



if __name__ == '__main__':
	collection = "../data/rome"
	collection_distributions_4in1(collection)
