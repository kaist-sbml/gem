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
RUN apt-get install -y sed

# Install major dependencies
COPY . /gmsm
WORKDIR /gmsm/
RUN pip install pip --upgrade
RUN pip install -r requirements.txt
RUN run_gmsm.py -p -s -d -i input/sample_input_two_CDS.gb
RUN rm -r /gmsm/output

# Set GMSM implementation
ADD ./EFICAz2.5.1.tar.gz /
WORKDIR /usr/bin/
RUN bash "/EFICAz2.5.1/bin/INSTALL"
RUN sed -i 's/\/usr/$EFICAz25_PATH/g' /usr/bin/eficaz2.5
VOLUME ["/input", "/output"]
WORKDIR /gmsm/
RUN rm /gmsm/EFICAz2.5.1.tar.gz
RUN chmod +x run_gmsm.py

ENTRYPOINT ["/bin/bash"]
