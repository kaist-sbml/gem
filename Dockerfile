# Install ubuntu
FROM ubuntu:14.04
MAINTAINER Hyun Uk Kim, Jae Yong Ryu

# Set environment variables
ENV EFICAz25_PATH="/EFICAz2.5.1/"
ENV PATH="${PATH}:${EFICAz25_PATH}"

ENV GMSM_PATH="/gmsm/"
ENV PATH="${PATH}:${GMSM_PATH}"

RUN apt-get -y update
RUN apt-get install -y python2.7
RUN apt-get install -y python2.7-dev
RUN apt-get install -y python-pip
RUN apt-get install -y python-tox
RUN apt-get install -y ncbi-blast+

# Install major dependencies
RUN pip install pip --upgrade
RUN pip install biopython==1.68
RUN pip install cobra>=0.6.2
RUN pip install pyparsing==2.2.0
RUN pip install pytest-cov==2.4.0
RUN pip install python-libsbml

# Set GMSM implementation
COPY . /gmsm
ADD ./EFICAz2.5.1.tar.gz /
WORKDIR /usr/bin/
RUN bash "/EFICAz2.5.1/bin/INSTALL"
VOLUME ["/input", "/output"]
WORKDIR /gmsm/
RUN chmod +x run_gmsm.py

ENTRYPOINT ["/bin/bash"]
