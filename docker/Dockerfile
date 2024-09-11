FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && apt-get install -y build-essential wget curl unzip cmake python3-pip mono-complete python-is-python3

COPY . /root/xic-extractor
WORKDIR /root/xic-extractor

RUN python -m pip install --user -r requirements.txt
ENV PATH="${PATH}://root/xic-extractor/bin"
