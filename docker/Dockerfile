FROM ubuntu:18.04

LABEL version="0.4.0" \
     description="APHA Bovine TB pipeline" \
     dockerhub="oxfordmmm/apha-btb:v0.4.0"

#Fresh update packages
RUN apt-get update && \
    apt-get install -y software-properties-common

RUN apt install -y \
       bc \
       git  \
       gcc \ 
       libncurses5-dev \
       libbz2-dev \
       liblzma-dev \
       make \
       python \
       python-numpy \
       python-pip \
       unzip \
       wget \
       zip \
       zlib1g-dev   &&\
    apt-get clean     &&\
    apt-get purge     &&\
    rm -rf /var/lib/apt/lists/ */tmp/* /var/tmp/*

# java8
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    apt install openjdk-8-jre -y

# biopython
RUN pip install biopython

# FastUniq
WORKDIR /biotools
RUN wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz &&\
    tar xzf FastUniq-1.1.tar.gz &&\
    rm -f FastUniq-1.1.tar.gz &&\
    mv FastUniq FastUniq-1.1 &&\
    cd FastUniq-1.1/source &&\
    make

# Trimmomatic
WORKDIR /biotools
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip &&\
    unzip Trimmomatic-0.38.zip &&\
    rm -f Trimmomatic-0.38.zip

# Bwa
WORKDIR /biotools
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 &&\
    tar -xvf bwa-0.7.17.tar.bz2 &&\
    rm bwa-0.7.17.tar.bz2 &&\
    cd bwa-0.7.17 &&\
    make 

# samtools
WORKDIR /biotools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 &&\
    tar xjf samtools-1.9.tar.bz2 &&\
    rm -f samtools-1.9.tar.bz2 &&\
    cd samtools-1.9 &&\ 
    make &&\
    make install

# bcftools
WORKDIR /biotools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 &&\
    tar xjf bcftools-1.9.tar.bz2 &&\
    rm -f bcftools-1.9.tar.bz2 &&\
    cd bcftools-1.9 &&\
    make &&\
    make install

# vcfutils.pl
WORKDIR /biotools
RUN wget https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl

#kraken2
WORKDIR /biotools
RUN wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz && tar xzf v2.0.8-beta.tar.gz && rm -f v2.0.8-beta.tar.gz
WORKDIR /biotools/kraken2-2.0.8-beta
RUN ./install_kraken2.sh ../Kraken2

#create btb user
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse aphabtb && \
    echo `echo "aphabtb\naphabtb\n" | passwd aphabtb` && \
    chown aphabtb:aphabtb /home/aphabtb

# script
WORKDIR /home/aphabtb
COPY pyscripts pyscripts

#env
ENV TOOLS=/biotools
ENV PATH=${PATH}:$TOOLS

# user
USER aphabtb
ENTRYPOINT [ "/bin/bash" ]

#To build this docker #docker build -t apha-btb:v0.4.0 .
#To run this docker   #docker run  -it apha-btb:v0.4.0
#To pull this docker  #docker pull oxfordmmm/apha-btb:v0.4.0
#To tranform to Singularity #singularity build apha-btb-v0.4.0.img docker://oxfordmmm/apha-btb:v0.4.0