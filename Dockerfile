FROM stackbrew/ubuntu:14.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git tar gzip bzip2 build-essential \
        python2.7-dev python-pip python-virtualenv zlib1g-dev
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
RUN apt-get install git-lfs
RUN git lfs install

# TargQC installation
RUN pip install --upgrade setuptools pip
RUN git clone --depth 50 https://github.com/vladsaveliev/TargQC.git TargQC
RUN cd TargQC && \
    git submodule update --init --recursive && \
    pip install -r requirements.txt && \
    python setup.py install


# clean filesystem
#    cd /usr/local && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/* /var/tmp/* && \
#    /usr/local/share/bcbio-nextgen/anaconda/bin/conda clean --yes --tarballs && \
#    rm -rf /usr/local/share/bcbio-nextgen/anaconda/pkgs/qt* && \
#    rm -rf /usr/local/.git && \
#    rm -rf /.cpanm && \
#    rm -rf /tmp/bcbio-nextgen-install && \

# Create directories and symlinks for data
#    mkdir -p /mnt/biodata && \
#    mkdir -p /tmp/bcbio-nextgen && \
#    mv /usr/local/share/bcbio-nextgen/galaxy/bcbio_system.yaml /usr/local/share/bcbio-nextgen/config && \
#    rmdir /usr/local/share/bcbio-nextgen/galaxy && \
#    ln -s /mnt/biodata/galaxy /usr/local/share/bcbio-nextgen/galaxy && \
#    ln -s /mnt/biodata/gemini_data /usr/local/share/bcbio-nextgen/gemini_data && \
#    ln -s /mnt/biodata/genomes /usr/local/share/bcbio-nextgen/genomes && \
#    ln -s /mnt/biodata/liftOver /usr/local/share/bcbio-nextgen/liftOver && \
#    chmod a+rwx /usr/local/share/bcbio-nextgen && \
#    chmod a+rwx /usr/local/share/bcbio-nextgen/config && \
#    chmod a+rwx /usr/local/share/bcbio-nextgen/config/*.yaml && \

# Ensure permissions are set for update in place by arbitrary users
#    find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
#    find /usr/local -perm /u+w -execdir chmod a+w {} \;
