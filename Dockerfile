FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update

RUN curl -o blast.tar.gz 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz'
RUN tar -zxvf blast.tar.gz

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/bin
RUN mv ./ncbi-blast-2.2.29+/bin/* /kb/module/bin/
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
