FROM kbase/sdkpython:3.8.10
MAINTAINER Dylan Chivian
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Fix Certs
RUN apt-get update
RUN apt-get upgrade -y
RUN sed -i 's/\(.*DST_Root_CA_X3.crt\)/!\1/' /etc/ca-certificates.conf
RUN update-ca-certificates 

RUN pip install --upgrade pip

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all


# MMseqs2
RUN curl --insecure -o mmseqs-linux-avx2.tar.gz https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    tar xvzf mmseqs-linux-avx2.tar.gz
ENV PATH="${PATH}:/kb/module/mmseqs/bin"

# mOTUlizer (needs upgraded pip - done above)
RUN pip install mOTUlizer


ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
