FROM ubuntu:15.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget make g++ git tar gzip bzip2 build-essential  \
        python2.7-dev python-pip python-virtualenv zlib1g-dev default-jre && \
    apt-get upgrade -y libstdc++6

# Install git LFS
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    apt-get install git-lfs && \
    git lfs install

# Matplotlib dependencies
RUN apt-get update && apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib

# TargQC installation
COPY . TargQC
RUN pip install --upgrade setuptools pip && \
    cd TargQC && \
    pip install --upgrade -r requirements.txt && \
    python setup.py develop && \
    cd ..
