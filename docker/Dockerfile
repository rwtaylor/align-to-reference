FROM debian:latest
MAINTAINER Ryan Taylor <ryan@ryantaylor.net>

RUN apt-get update -y && apt-get install -y apt-utils && apt-get dist-upgrade -y 
RUN apt-get install -y build-essential zlib1g-dev wget default-jre unzip git python python-pip \
    python-dev libncurses5-dev libgsl0-dev pkg-config libfreetype6-dev libpng-dev && \
    apt-get autoremove

# seqtk
#RUN apt-get install -y  zlib1g-dev
RUN cd /usr/local/src && \
    git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && git checkout tags/v1.2 && \
    make && cp seqtk /usr/local/bin/

#fastqc / fastqt
RUN mkdir -p /usr/local/opt && \
    cd /usr/local/opt && \
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    cd FastQC && chmod +x fastqc && \
    ln -s /usr/local/opt/FastQC/fastqc /usr/local/bin/fastqc

# cutadapt
#RUN apt-get install -y python python-pip python-dev 
RUN pip install cutadapt

#trim_galore
RUN mkdir -p /usr/local/src/trim_galore && \
    cd /usr/local/src/trim_galore && \
    wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.3.zip && \
    unzip trim_galore_v0.4.3.zip && \
    cp trim_galore /usr/local/bin/trim_galore

#bwa mem
RUN cd /usr/local/src && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa && git checkout tags/v0.7.15 && \
    make && \
    cp bwa /usr/local/bin/bwa

#samtools
#RUN apt-get install -y libncurses5-dev
RUN cd /usr/local/src && \
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xvf samtools-1.* && \
    cd samtools-1.* && \
    make && \
    make prefix=/usr/local install

#picardtools
RUN cd /usr/local/opt && \
    wget https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar

# htslib
RUN cd /usr/local/src && \
    wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar -xvf htslib-1.3.2.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && make install

#preseq
#RUN apt-get install -y libgsl0-dev
RUN cd /usr/local/src && \
    wget http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 && \
    tar -xvf preseq_linux_v2.0.tar.bz2 && \
    cd preseq_v2.0 && \
    mv preseq /usr/local/bin/ && mv bam2mr /usr/local/bin/
    
#multiqc
#RUN apt-get install -y pkg-config libfreetype6-dev libpng-dev
RUN pip install multiqc


RUN echo "deb http://ftp.debian.org/debian testing main" >> /etc/apt/sources.list && \
    apt-get update -y && apt-get install -y openjdk-8-jre