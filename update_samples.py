import yaml
config_path = "/Volumes/SSD/plan_a/visiumHD_pipeline_2/config/pipeline.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Update the condition test default samples to 50
if "condition_test" not in config:
    config["condition_test"] = {}
config["condition_test"]["samples"] = 50
config["condition_test"]["recorded_samples"] = 20

with open(config_path, "w") as f:
    yaml.dump(config, f, default_flow_style=False)
