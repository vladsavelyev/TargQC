FROM stackbrew/ubuntu:14.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git tar gzip bzip2 build-essential \
        python2.7-dev python-pip python-virtualenv zlib1g-dev
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    apt-get install git-lfs && \
    git lfs install

# TargQC installation
RUN pip install --upgrade setuptools pip
RUN git clone --recursive https://github.com/vladsaveliev/TargQC.git TargQC
RUN cd TargQC && \
    git submodule update --init --recursive && \
    pip install --upgrade -r requirements.txt && \
    python setup.py develop
