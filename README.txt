This is the Snakemake Pipeline for (1) simulating mutispecies predator prey dynamics with gLV and (2) inferring networks with 7 inference Methods and (3) Analysis script to assess inference qualities. This requires working Python and R environments, code Editors like Spyder and Rstudio are not strictly required but recommended to make changes to the scripts. Here, the Setup process of R requirements is described using Rstudio.

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

install.packages("NetComi")

for all listed packages.
-------

