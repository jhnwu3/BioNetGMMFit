# os to use for building
FROM ubuntu AS build
# build chain that needs to happen
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential git cmake autoconf libtool pkg-config 

COPY /mnt/c/buildroadrunner /buildroadrunner
COPY . /CyGMM/
WORKDIR "/CyGMM"
RUN cd src && \
	cmake . && \
	make && \
	cd .. 

# os to use for running, which is still big.
FROM ubuntu
COPY --from=build /CyGMM /CyGMM
WORKDIR "/CyGMM"
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y python3 python3-pip && \
    pip install bionetgen
	
ENTRYPOINT [ "./CyGMM"]