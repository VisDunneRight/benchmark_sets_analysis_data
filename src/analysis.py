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


def collection_distributions(path_to_folder, create_dir):
	fpath = path_to_folder + "/clean"
	node_counts = []
	edge_counts = []
	avg_degree_counts = []
	max_degree_counts = []
	for gfile in os.listdir(fpath):
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
	if create_dir:
		os.mkdir(path_to_folder + "/charts")
	plot_histogram(node_counts, "Node Count", path_to_folder + "/charts/node_counts.svg")
	plot_histogram(edge_counts, "Edge Count", path_to_folder + "/charts/edge_counts.svg")
	plot_histogram(avg_degree_counts, "Mean Degree", path_to_folder + "/charts/average_degree.svg")
	plot_histogram(max_degree_counts, "Maximum Degree", path_to_folder + "/charts/max_degree.svg")


if __name__ == '__main__':
	collection_distributions("../data/rome", False)
