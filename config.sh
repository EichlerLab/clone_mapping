module purge

# Setup Python environment.
export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis
export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts
export PYTHONPATH=$PYTHONPATH:/net/gs/vol2/home/psudmant/EEE_Lab/projects/common_code
export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/projects/common_code/ssf_DTS_caller
export PYTHONPATH=$PYTHONPATH:/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/get_gc_correction

module load modules modules-init modules-gs/prod modules-eichler/prod

module load samtools/1.3 bedtools/2.23.0 mrsfast/3.3.8 anaconda/2.3.0 ucsc/20140617 openmpi/1.5.4 tabix/0.2.6

