FROM ubuntu:14.04

MAINTAINER Travis Hoppe

ENV JAVA_VERSION_MAJOR 8
ENV JAVA_VERSION_MINOR 66
ENV JAVA_VERSION_BUILD 17
ENV JAVA_PACKAGE       server-jre

ADD matlab.txt /mcr-install/matlab.txt

RUN apt-get update && \
	apt-get install -y curl wget unzip xorg

# Download and unarchive Oracle Java
RUN curl -jksSLH "Cookie: oraclelicense=accept-securebackup-cookie"\
  http://download.oracle.com/otn-pub/java/jdk/${JAVA_VERSION_MAJOR}u${JAVA_VERSION_MINOR}-b${JAVA_VERSION_BUILD}/${JAVA_PACKAGE}-${JAVA_VERSION_MAJOR}u${JAVA_VERSION_MINOR}-linux-x64.tar.gz \
    | tar -xzf - -C /opt &&\
    ln -s /opt/jdk1.${JAVA_VERSION_MAJOR}.0_${JAVA_VERSION_MINOR} /opt/jdk

# Install MatLab runtime
RUN \
	cd /mcr-install && \
	wget -nv http://www.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip && \
	unzip MCR_R2013a_glnxa64_installer.zip && \
    cd /mcr-install && \
	mkdir /opt/mcr && \
	./install -inputFile matlab.txt && \
	cd / && \
	rm -rf mcr-install

# Set environment
ENV JAVA_HOME /opt/jdk
ENV PATH ${PATH}:${JAVA_HOME}/bin

# ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/mcr/v81/bin/glnxa64/:/opt/mcr/v81/runtime/glnxa64

ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/mcr/v81/runtime/glnxa64:/opt/mcr/v81/bin/glnxa64:/opt/mcr/v81/sys/os/glnxa64:/opt/mcr/v81/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/opt/mcr/v81/sys/java/jre/glnxa64/jre/lib/amd64/server:/opt/mcr/v81/sys/java/jre/glnxa64/jre/lib/amd64



# Define default command
CMD ["bash"]


