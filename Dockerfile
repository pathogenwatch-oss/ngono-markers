FROM python:3.10-slim

ARG DEBIAN_FRONTEND=noninteractive

COPY --from=docker.io/astral/uv:latest /uv /uvx /bin/

RUN apt update \
    && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    build-essential \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /data

COPY db/markers.fasta /db/

RUN cd /db && makeblastdb -in markers.fasta -out markers -dbtype nucl -parse_seqids

COPY marker_search /marker_search

COPY uv.lock pyproject.toml ngono-markers.py /

RUN uv sync --locked && \
    uv run /ngono-markers.py --help

#CMD ["/bin/bash"]
ENTRYPOINT ["uv", "run", "/ngono-markers.py", "--blast-db", "/db"]