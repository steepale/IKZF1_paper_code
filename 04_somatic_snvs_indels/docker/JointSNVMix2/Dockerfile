#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it --name jsm2 steepale/jsm2:0.7.5
#        USAGE: docker image build -t steepale/jsm2:4.0.2 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build JointSNVMix2
# REQUIREMENTS:  ---
#        NOTES:  Instructions: https://www.gungorbudak.com/blog/2014/03/26/jointsnvmix-installation-on-linux-mint/ (pysam 0.5.0)
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2019.09.26
#     REVISION:  ---
#===============================================================================

# Pull the Ubuntu OS image
FROM steepale/ubuntu_base:18.04

# Change to base directory
WORKDIR /

# Install the jsm2 dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN sudo apt update \
    && sudo apt -y dist-upgrade \
    && sudo apt install -y python2.7 \
    python-dev \
    libevent-dev \
    python-pip \
    zlib1g-dev \
    && sudo pip install numpy

# Install application software in opt directory
WORKDIR /opt

# Install Cython 0.20.1
RUN wget https://files.pythonhosted.org/packages/81/87/9ceffc2c15a06fcdd82e621b54598da684271ed0c6722b316e7a30e4c18e/Cython-0.20.1.tar.gz \
    && tar zxvf Cython-0.20.1.tar.gz \
    && cd Cython-0.20.1 \
    && sudo python setup.py install


# Install Anaconda (requires interaction)
#RUN wgte https://repo.anaconda.com/archive/Anaconda2-2019.07-Linux-x86_64.sh \
#    && sudo bash Anaconda2-2019.07-Linux-x86_64.sh \
#    && echo 'export PATH="/root/anaconda2/bin:${PATH}"' >> ~/.bashrc

# Install Cython (http://www.alglib.net/download.php)
RUN wget http://www.alglib.net/translator/re/alglib-3.15.0.cpython.free.tgz \
    && tar zxvf alglib-3.15.0.cpython.free.tgz \
    && rm alglib-3.15.0.cpython.free.tgz \
    && cd cpython \
    && python setup.py install

# Install application software in opt directory
WORKDIR /opt

# To install pysam, first install distribute
RUN wget https://pypi.python.org/packages/source/d/distribute/distribute-0.7.3.zip \
    && unzip distribute-0.7.3.zip \
    && rm distribute-0.7.3.zip \
    && cd distribute-0.7.3 \
    && sudo python setup.py install \
    && cd ..

# Then install ez_setup
RUN wget https://pypi.python.org/packages/source/e/ez_setup/ez_setup-0.9.tar.gz \
    && tar xzvf ez_setup-0.9.tar.gz \
    && rm ez_setup-0.9.tar.gz \
    && cd ez_setup-0.9 \
    && sudo python setup.py install \
    && cd ..

# Adjust dependencies in Docker file

sudo apt-get install -y gcc-4.8
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50


# Install pysam 0.5.0 (keep this tar file, it's no longer available)
RUN wget http://pysam.googlecode.com/files/pysam-0.5.tar.gz \
    && tar xzvf pysam-0.5.tar.gz \
    && cd pysam-0.5 \
    && sudo python setup.py install \
    && cd ..

# Install jsm2
RUN wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/joint-snv-mix/JointSNVMix-0.7.5.tar.gz \
    && tar zxvf JointSNVMix-0.7.5.tar.gz \
    && rm JointSNVMix-0.7.5.tar.gz \
    && cd JointSNVMix-0.7.5 \
    && python setup.py install \
    && echo 'export PATH="/opt/JointSNVMix-0.7.5/include/pysam/:${PATH}"' >> ~/.bashrc \


# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]