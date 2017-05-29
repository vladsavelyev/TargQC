# Builds and uploads packages for macos and linux, py2 and py3

function build() {
	local NAME=$1
	CHANNELS="-c vladsaveliev -c bcbio -c bioconda -c conda-forge"
	DIR=conda-bld
	mkdir -p $DIR
	PATH_36=$(conda build $NAME $CHANNELS --output-folder $DIR --output --python 3.6 | tail -n1)
	PATH_27=$(conda build $NAME $CHANNELS --output-folder $DIR --output --python 2.7 | tail -n1)
	conda convert -p linux-32 $PATH_36
	conda convert -p linux-64 $PATH_36
	conda convert -p linux-32 $PATH_27
	conda convert -p linux-64 $PATH_27
	FILE_36=$(basename $PATH_36)
	FILE_27=$(basename $PATH_27)
	for p in osx-64 linux-32 linux-64
	do
		anaconda upload $p/$FILE_36
		anaconda upload $p/$FILE_27
	done
}

# Iterate over directories and build all packages
for f in */meta.yaml
do 
	package_name=$(dirname ${f#/*});
	echo "Building $package_name";
	build $package_name;
done
