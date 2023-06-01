import json
import os


def read_json(path_to_file):
	with open(path_to_file) as fp:
		gr_json = json.load(fp)
	return gr_json


def collection_distributions(path_to_folder):
	fpath = path_to_folder + "/clean"
	node_counts = []
	edge_counts = []
	for gfile in os.listdir(fpath):
		with open(fpath + "/" + gfile) as fd:
			gr = json.load(fd)
		node_counts.append(len(gr["nodes"]))
		edge_counts.append(len(gr["links"]))
