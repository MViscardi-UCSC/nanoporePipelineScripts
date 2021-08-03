#FROM ubuntu:xenial
FROM nvidia/cuda:11.1-cudnn8-devel-ubuntu16.04
USER 0
RUN apt update -y
RUN apt upgrade -y

COPY python3-ont-albacore_2.3.4-1_xenial_amd64.deb albacore.deb
COPY ont_guppy_cpu_4.2.2-1~xenial_amd64.deb guppy.deb
# libboost is important for C+ functionality (I think?)
RUN apt install libboost-all-dev apt-utils -y
ENV DEBIAN_FRONTEND=noninteractive
# This installs 3.5 with xenial, but 3.7 in bionic
RUN apt install python3-pip -y
RUN pip3 install --upgrade pip

# Albacore Needs + Install:
###########################
RUN apt install libcurl3 libzmq5 libhdf5-dev -y
# Trying to setup a fake python3-ont-fast5-api so that albacore realizes its being ridiculous and can stop saying its not installed
COPY python3-ont-fast5-api_3.3_all.deb fast5.deb
RUN dpkg -i fast5.deb
RUN python3 -m pip install ont-fast5-api mappy
# Finish Albacore stuff
RUN apt install python3-h5py python3-numpy  python3-dateutil python3-progressbar python3-setuptools -y
RUN dpkg --force-confold -i albacore.deb

# Guppy Install
###############
RUN apt install libcurl4-openssl-dev libssl-dev -y 
RUN dpkg --force-confold -i guppy.deb

# Nanopolish Install + Minimap2 + and maybe others...
#####################################################
RUN apt install git wget -y
RUN git clone --recursive https://github.com/jts/nanopolish.git
RUN cd nanopolish && make
RUN ln /nanopolish/nanopolish /usr/bin
RUN cd nanopolish/minimap2 && make
RUN ln /nanopolish/minimap2/minimap2 /usr/bin
RUN apt install samtools

# Bash starts up when using 'docker run -it <name>'
RUN apt install nano
RUN mkdir /usr/src/working
WORKDIR /usr/src/working

# Setup for specifically working with Roach et al. data:
# NOTE: This is set up to specifically run within /data16/marcus/working/210106_nanopolish_wRoachData  -  DO NOT MOVE IT!
RUN apt install tree
COPY docker_script.sh /bin/startup
CMD /bin/bash startup
