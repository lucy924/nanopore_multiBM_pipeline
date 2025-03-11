from jinja2 import Template
import codecs
import snakemake
import json

# This was from https://stackoverflow.com/a/63622717

# ------------------------------------------------ #
# # get snakemake variables
# sample_name = snakemake.params['sample_name']
# meth_results_fp = snakemake.input["epic_probe_results"]
# deconv_fp = snakemake.input["methatlas_deconvolution"]
# results_json_imin = snakemake.output["results_json_imin"]

# # path2_panel_metadata_csv = snakemake.input["panel_metadata"]
# # results_fp = snakemake.output['panel_mod_results']

# log = open(snakemake.log[0], 'w')

# ------------------------------------------------ #
sample_name = "test"

# create an dict will all data that will be populate the template
path2results_dict = (
    "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/results/test2/results.imin.json"
)

with open(path2results_dict, "r") as fj:
    results_dict = json.load(fj)

results = {
    "name": sample_name,
    "url": "https://github.com/lucy924/Multi-biomarker-ONT-project",
}

results["results"] = list_of_lists = [
    [key, value] for key, value in results_dict.items()
]

# render the template
path2template = (
    "/home/dejlu879/ProjectProtocol/ext_bm_pipeline_dev/resources/template.md"
)
with open(path2template, "r") as file:
    template = Template(file.read(), trim_blocks=True)
rendered_file = template.render(repo=results)

# output the file
output_file = codecs.open("report.md", "w", "utf-8")
output_file.write(rendered_file)
output_file.close()
