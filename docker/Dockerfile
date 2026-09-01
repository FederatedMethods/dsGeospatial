FROM datashield/rock-base:6.3.0-R4.4.0

COPY dsGeospatial_0.0.0.9000.tar.gz /tmp/

RUN R CMD INSTALL /tmp/dsGeospatial_0.0.0.9000.tar.gz

RUN rm /tmp/dsGeospatial_0.0.0.9000.tar.gz
