This is the Snakemake Pipeline for (1) simulating mutispecies predator prey dynamics with gLV and (2) inferring networks with 7 inference Methods and (3) Analysis script to assess inference qualities. This requires working Python and R environments, code Editors like Spyder and Rstudio are not strictly required but recommended to make changes to the scripts. Here, the Setup process of R requirements is described using Rstudio.

Workflow description:

generateBasicSyntheticNetwork - generates a synthetic network topology (cluster, scale free, band) with n-species, which serves as the basis for simulations
generateSyntheticAbundaceData_gLV - edges of synthetic network get assigned interactions randomly, then the species abundances are simulated with the generalized Lotka-Volterra model. Attractor states serve as a sample, contributing to the final abundance Matrix. In this step ESABO network is inferred from the abundances.
inferNetworksR - given the simulated abundances, networks are inferred using the R environment (SpiecEasi, CCREPE, SPARCC, Spearman, propr, ecocopula) 
benchmarkInferenceQuality - inference results are aggreated and qualitatively compared to the base synthetic network. Plots Show the TP/FP rate and 

++ 1. Setup Python Environment 

Python & Package Management e.g. Miniforge:
https://github.com/conda-forge/miniforge
make sure it works in the terminal:
conda --Version

Clone this repository into a Project Folder:
cd /your/project/folder
git clone https://github.com/v-dinkel/FoodWeb_gLV

Change to the new subfolder /your/project/Folder/FoodWeb_gLV
Create new Environment from given Project requirements. Other versions may also work but this is the verified Environment:

conda create --name foodweb_glv --file requirements.txt
conda activate foodweb_glv

if conda is not working, try this command:
source ~/.bashrc

The requirements doesn't contain snakemake, so install it as stated here:
https://gist.github.com/RomainFeron/da9df092656dd799885b612fedc9eccd

or use this command:
conda install -c conda-forge -c bioconda snakemake

make sure ist working with snakemake --Version

Using an editor like Spyder, set the Interpreter to the Environment. Go to 
Tools > Preferences > Python Interpreter > Use the following Interpreter > /home/user/miniforge3/envs/foodweb_glv/bin/python

++ 2. R Environment

Install R and RStudio. If you are using Ubuntu.

Start RStudio and open the script inferNetworks.R. R Studio should automatically suggest the Installation of missing packages. If this does not work, you can manually install the packages using:

install.packages("devtools")
install.packages("BiocManager")
...

Some packages (devtools) require System dependencies. This is highly individual and system dependant, therefore see which 'ERROR: dependency '... is missing and try solving it. Also, the packages SpiecEasi, SPRING and NetComi require an extensive list of packages, which also may be in Conflict with version of the OS and R. Package version dependencies are best adressed with:

install.packages("remotes")
remotes::install_version("Matrix", version="1.6")

++ 3. Configure Snakemake
Change the workdir in the config.yaml to your folder. This path must be absolute.
Other variables are:
n is the amount of species to be simulated. This generates a nxn matrix. Experiments were performed with n=20. Large n significantly increases processing speed.
seed is the random state.
nettype defines the topology of the network which will be simulated. They can be cluster, scale_free and band. This uses the SpiecEasi package, for more options please see SpiecEasi documentation.
r_methods are the networks which need to be inferred. 6 methods are in R, ESABO is inferred during the simulation of abundance data.
nsimulations is the amount of abundance simulations and corresponding network inferrences.

## 4. Run Snakemake Pipeline

The Pipeline doesn't use snakemake rule dependencies to accomplish a fully automated workflow. This is intentional since simulations and network inferences are time consuming and large benchmarks can have unexpected Errors at some Point in the execution. This is to avoid snakemake cleanup process and keep all generated files for error inspection or partial analyses.

Run the individual steps by their Name. The -j 4 command defines the amount of processors used. Adapt to your Settings.

snakemake -j 4 generateBasicSyntheticNetwork
snakemake -j 4 generateSyntheticAbundaceData_gLV
snakemake -j 4 inferNetworksR
snakemake -j 4 benchmarkInferenceQuality

## 5. Outputs

outputs/<seed>/

abundances - contains abundance files from n-simulations. Plots illustrate an example sample Simulation (generalized lotka volterra) where only the last values (attractors) are extracted. An abundance table is formed by attractors (
benchmark - aggregated results of the simulations and plots
Graphs - synthetic topology graphs
networks - inferred networks for each simulation


Notes:

- Make sure to remove/comment package installation commands from the R script
- In the first execution, some missing R packages may be installed by NetComi.
- During network inference, sometimes the text output in the terminal appears to be stuck (no new lines). If that happens, try pressing Enter in the terminal.
