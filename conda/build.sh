# Builds and uploads packages for macos and linux, py2 and py3

set -x

function build() {
	local NAME=$1
	CHANNELS="-c vladsaveliev -c bcbio -c bioconda -c conda-forge"

	for PY in 3.6 2.7 
	do
		PACKAGE_PATH=$(conda build $NAME $CHANNELS --output --py $PY | tail -n1)
		echo "Building $PACKAGE_PATH"
		conda build $NAME $CHANNELS --py $PY
		BASEDIR=$(dirname $(dirname $PACKAGE_PATH))		
		FILENAME=$(basename $PACKAGE_PATH)
		echo "Converting packages into $BASEDIR"
		for PLATFORM in osx-64 linux-32 linux-64
		do
			conda convert -p $PLATFORM $PACKAGE_PATH -o $BASEDIR
			anaconda upload $BASEDIR/$PLATFORM/$FILENAME
		done
	done
	cd ..
}

# Iterate over directories and build all packages
for f in */meta.yaml
do 
	package_name=$(dirname ${f#/*});
	echo "Building $package_name";
	build $package_name;
done

set +x