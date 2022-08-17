FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:6839-main

RUN apt-get install -y curl unzip

RUN apt install -y gfortran
RUN python3 -m pip install scanpy numpy pandas 
RUN python3 -m pip install dropkick matplotlib leidenalg nbformat jupyter nbconvert

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch

WORKDIR /root
