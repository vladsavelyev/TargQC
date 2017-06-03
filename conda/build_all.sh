#!/usr/bin/env bash
# Builds and uploads packages for linux-32, linux-64, and macos, for py2 and py3

set -x
set -e

function build() {
	local NAME=$1
	echo "Building $NAME";

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
}

if [ -z "$1" ]; then
	# Iterate over directories and build all packages
	for f in */meta.yaml
	do
		package_name=$(dirname ${f#/*});
		build $package_name;
	done
else
	build $1;
fi


set +x
set +e
