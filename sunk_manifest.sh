#!/bin/env bash

while getopts ":d:o:s:" opt; do
	case ${opt} in
		d )
		directory=$OPTARG
		;;
		o )
		outfile=$OPTARG
		;;
		s )
		samplefile=$OPTARG
		;;
		\? )
		echo "Usage: ./sunk_manifest.sh [-d {directory where to look for files} (default is ./)] 
		[-s SampleSheet.csv containing pool info (usually contained 3 directories up from the output fastq files directory)]  
		[-o {output file name}]"
		exit 0
		;;
	esac
done

if [ -z ${directory+x} ] 
	then 
		echo "directory is unset please specify with -d"
		exit 1
	else 
		echo -e "Checking ${directory} for SUNK manifest"
fi

if [ -z ${outfile+x} ] 
	then 
		echo "Output file is unset please specify with -o"
		exit 1
	else 
		echo -e "Output will be directed to ${outfile}"
fi

if [ -z ${samplefile+x} ] 
	then 
		echo "SampleSheet.csv is unset please specify with -s"
		exit 1
	else 
		echo -e "SampleSheet.csv found"
fi


tmp_directory=${RANDOM}.dir

/bin/ls ${directory}/*fastq | sort > ${tmp_directory}

clone_names=$( cat ${tmp_directory} | awk -F "/" '{print $NF}' | awk -F '_' '{print $1}' | sort -u )

echo -e "well\tsample_name\treads" > ${outfile}

for fn in ${clone_names}; do 
	sn=$( grep ${fn} ${tmp_directory} | head -1 ) 
	sn2=$( grep ${fn} ${tmp_directory} | tail -1 )
	well=$( grep ${fn} ${samplefile} | awk -F "," '{print $NF}' )
	echo -e "${well}\t${fn}\t\"${sn},${sn2}\"" | sed -e "s/\r//g" >> ${outfile}.tmp
done

sort -V ${outfile}.tmp >> ${outfile}

rm ${outfile}.tmp
rm ${tmp_directory}
