# os to use
FROM ubuntu:20.04 

# RUN chain to install all dependencies (?) first let's see if we can't just run it with whatever we have now without linked libraries
# RUN apt-get update && \
#     apt-get install gcc

RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config
    
COPY . /CyGMM 
ENTRYPOINT [ "CyGMM/CyGMM" ]
