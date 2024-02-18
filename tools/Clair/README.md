### Installation

```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "clair-env"
conda create -n clair-env -c bioconda clair
conda activate clair-env
conda install deepdish

# store clair.py PATH into $CLAIR variable
CLAIR=`which clair.py`

# run clair like this afterwards
python $CLAIR --help
```

The conda environment has the Pypy3 interpreter installed, but one Pypy3 package `intervaltree` is still missing. The reason why this is not installed by default is because this is not yet available in any conda repositories. To install the package for Pypy3, after activating the conda environment, please run the following commands:

```bash
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
```

### Create dataset

```shell
python3 clair.py callVarBam --chkpnt_fn ./model --bam_fn <chr.bam.bai> --ref_fn <chr20.fa>  --call_fn .chr20.vcf --pysam_for_all_indel_bases --threads 1 --qual 100 --ctgName chr20 --ctgStart 10269870 --ctgEnd 46672937 --pipe_line --store_loaded_mini_match
```

### Run command

Run command
```shell
# download models
wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
tar -xf 12345.tar
python3 my_prediction.py --chkpnt_fn ./model --sampleName HG001 --threads 1 --qual 100 --input_fn <prediction_input.h5> --output_fn <prediction_output.h5>
```
