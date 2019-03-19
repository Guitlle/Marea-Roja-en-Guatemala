# docker build -t geopy-1 . 

FROM python:3.7-stretch

RUN apt-get update
RUN apt-get install libgeos-3.5.1 libgeos-c1v5 libgeos-dev

COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

RUN mkdir /rawdata
WORKDIR /
#COPY producto /app

# to test the container in bash:
#    docker run -it --rm geopy-1 bash
# to run the fan script:
#    docker run -it --rm -v `pwd`/producto:/app -v `pwd`/output:/output geopy-1 python3 /app/fan_gt.py
