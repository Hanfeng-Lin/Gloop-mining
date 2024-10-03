FROM ubuntu:18.04

# Install dependencies
RUN apt-get update
RUN apt-get install -y wget python3 python3-pip python3-venv vim dssp
RUN pip3 install --upgrade pip && \
    pip3 install scipy biopython