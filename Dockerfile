FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
    && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install biopython click

RUN mkdir /data

COPY db/markers.fasta /db/

RUN cd /db && makeblastdb -in markers.fasta -out markers -dbtype nucl -parse_seqids

COPY marker_search /marker_search

COPY ngono-markers.py /

#CMD ["/bin/bash"]
ENTRYPOINT ["python3", "/ngono-markers.py"]