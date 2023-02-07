# SOURCE AND CREDITS: https://github.com/chrishah/samtools-docker/blob/master/Dockerfile

FROM alpine:latest
LABEL version="v1.0" description="A docekr image for the Pharmacogenetics ANalysis Pipeline" maintainer="Graeme Ford <Graeme.Ford@tuks.co.za>"

RUN apk update && apk upgrade

# Install python/pip
ENV PYTHONUNBUFFERED=1
RUN apk add --update --no-cache python3 && ln -sf python3 /usr/bin/python
RUN python3 -m ensurepip
RUN pip3 install --no-cache --upgrade pip setuptools

RUN apk add --no-cache  build-base 
RUN apk add --no-cache  ncurses-dev 
RUN apk add --no-cache  zlib-dev 
RUN apk add --no-cache  bzip2-dev
RUN apk add --no-cache  xz-dev

RUN apk add --no-cache  curl

WORKDIR /tmp

#Samtools
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools
ENV SAMTOOLS_VERSION=1.16.1

RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    tar --bzip2 -xf samtools-$SAMTOOLS_VERSION.tar.bz2
WORKDIR /tmp/samtools-$SAMTOOLS_VERSION
RUN ./configure --prefix=/usr && \
    make && \
    make install

# BCFTools
WORKDIR /tmp
ENV BCFTOOLS_VERSION=1.16

RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
    tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2
WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN ./configure -prefix=/usr/ && \
    make && \
    make install

# HTSLib
WORKDIR /tmp
ENV HTSLIB_VERSION=1.16

RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
    tar --bzip2 -xf htslib-$HTSLIB_VERSION.tar.bz2
WORKDIR /tmp/htslib-$HTSLIB_VERSION
RUN ./configure -prefix=/usr/ && \
    make && \
    make install

WORKDIR /work
# Plink-2.0
RUN wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip
RUN unzip plink2_linux_x86_64_20221024.zip
RUN mv plink2 /usr/bin/

ENV PATH="${PATH}:/usr/bin/"
