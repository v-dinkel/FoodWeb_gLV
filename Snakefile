configfile: "config.yaml"

rule generateBasicSyntheticNetwork:
	priority: 100
	output:
		"outputs/"+str(config['seed'])+"/graphs/base_graph_"+config['nettype']+".csv",
	script:
		"scripts/1_generateBaseNetworks.R"

rule generateSyntheticAbundanceData_gLV:
	priority: 99
	input:
		"outputs/"+str(config['seed'])+"/graphs/base_graph_"+config['nettype']+".csv",
	output:
		expand("outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}.csv", sim=[str(k) for k in range(config['nsimulations'])]),
		expand("outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_abunds.csv", sim=[str(k) for k in range(config['nsimulations'])]),
		expand("outputs/"+str(config['seed'])+"/networks/glv_"+config['nettype']+"_{sim}_esabo.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	script:
		"scripts/2_glv_simulation.py"

rule inferNetworksR:
	priority: 98
	input:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_abunds.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	output:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/networks/glv_"+config['nettype']+"_{sim}_{r_method}.csv", sim=[str(k) for k in range(config['nsimulations'])], r_method=[k for k in config['r_methods']]),
	script:
		"scripts/3_inferNetworks.R"
		
rule benchmarkInferenceQuality:
	priority: 97
	input:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_base_A.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	output:
		config['workdir']+"outputs/"+str(config['seed'])+"/benchmark/method_synthetic_benchmark.csv",
		config['workdir']+"outputs/"+str(config['seed'])+"/benchmark/benchmark_spieceasi_initial_cn.png"
	script:
		"scripts/4_qualityAssessment.py"

rule inferOnlySENetworksR:
	priority: 96
	input:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_abunds.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	output:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/networks/glv_"+config['nettype']+"_{sim}_spieceasi_weighted.csv", sim=[str(k) for k in range(config['nsimulations'])], r_method=['spieceasi']),
	script:
		"scripts/5_inferNetworks_SE_HT.R"