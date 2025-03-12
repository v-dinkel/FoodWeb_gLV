configfile: "config.yaml"

'''
rule benchmarkResults:
	input:
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free_partyHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_partyHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_dateHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band_partyHub.csv"
	output:
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free_partyHub.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_partyHub.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_dateHub.js",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band_partyHub.js"
	script:
		"scripts/benchmarkResults.py"

rule inferNetworksR:
	input:
		"outputs/synthetic_data/abundances/abund_scale_free.csv",
		"outputs/synthetic_data/abundances/abund_cluster.csv",
		"outputs/synthetic_data/abundances/abund_band.csv",
		"outputs/synthetic_data/abundances/abund_scale_free_partyHub.csv",
		"outputs/synthetic_data/abundances/abund_cluster_partyHub.csv",
		"outputs/synthetic_data/abundances/abund_cluster_dateHub.csv",
		"outputs/synthetic_data/abundances/abund_band_partyHub.csv"
	output:
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_scale_free_partyHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_partyHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_cluster_dateHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopula_band_partyHub.csv",
		"outputs/synthetic_data/inferred_graphs/ecoCopulaAUCs.csv"
	script:
		"scripts/inferNetworks.R"

rule DgenerateSyntheticAbundanceData_KL77:
	input:
		expand("outputs/synthetic_data/graphs/graph_{nettype}.csv", nettype=[k for k in config['nettypes']]),
	output:
		expand("outputs/synthetic_data/abundances/abund_kl77_{rank}_{nettype}.csv", rank=[k for k in config['ranks']], nettype=[k for k in config['nettypes']]),
	script:
		"scripts/generateSyntheticAbundance_KL77.R"

rule CgenerateSyntheticAbundanceData_Gut:
	input:
		expand("outputs/synthetic_data/graphs/graph_{nettype}.csv", nettype=[k for k in config['nettypes']]),
	output:
		expand("outputs/synthetic_data/abundances/abund_gut_{nettype}.csv", nettype=[k for k in config['nettypes']]),
	script:
		"scripts/generateSyntheticAbundance_Gut.R"
'''
rule generateBasicSyntheticNetwork:
	priority: 100
	output:
		"outputs/"+str(config['seed'])+"/graphs/base_graph_"+config['nettype']+".csv",
	script:
		"scripts/1_generateBaseNetworks.R"

rule GenerateSyntheticAbundanceData_gLV:
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
		
rule assessInferenceQuality:
	priority: 97
	input:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_base_A.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	script:
		"scripts/4_qualityAssessment.py"

rule inferOnlySENetworksR:
	priority: 96
	input:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/abundances/glv_"+config['nettype']+"_{sim}_filt_abunds.csv", sim=[str(k) for k in range(config['nsimulations'])]),
	output:
		expand(config['workdir']+"outputs/"+str(config['seed'])+"/networks/glv_"+config['nettype']+"_{sim}_{r_method}_HT.csv", sim=[str(k) for k in range(config['nsimulations'])], r_method=['spieceasi']),
	script:
		"scripts/5_inferNetworks_SE_HT.R"





