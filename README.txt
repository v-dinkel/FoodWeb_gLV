Install Python Environment & Package Management e.g. Miniforge:
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

The requirements doesn't contain snakemake, so install it as stated here:
https://gist.github.com/RomainFeron/da9df092656dd799885b612fedc9eccd

or use this command:
conda install -c conda-forge -c bioconda snakemake
-------

