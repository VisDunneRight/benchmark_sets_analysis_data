import json


def read_json(path_to_file):
	with open(path_to_file) as fp:
		gr_json = json.load(fp)
	return gr_json
