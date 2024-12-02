FROM rocker/r-ver:4.2.3
LABEL org.opencontainers.image.authors="viktor.s.petuhov@ya.ru"

# RUN apt-get  update --yes && apt-get install -y gnupg

# RUN apt-get remove libappstream3
# RUN apt-get --allow-insecure-repositories update --yes 
# RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 871920D1991BC93C

RUN apt-get && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libz-dev \
  nano \
  pandoc

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools", "Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown", "withr"), dependencies = TRUE)'

RUN mkdir /project
WORKDIR "/project"

RUN \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/4.2.3

# RUN \
#   R -e 'chooseCRANmirror(ind=52); install.packages(c("rlang", "httr", "gargle", "vctrs", "devtools"), dependencies = TRUE)'

RUN \
  R -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=T)' && \
  R -e 'chooseCRANmirror(ind=52); install.packages("ks", dependencies=c("Depends", "Imports", "LinkingTo"))' && \
  R -e 'devtools::install_local("~/dropEst/dropestr/", dependencies=T)'

RUN \
  mkdir -p dropEst/build && \
  cd dropEst/build && \
  cmake ../ && \
  make

ENV PATH=/home/user/dropEst/build:$PATH \
  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH \
  R_PROFILE=~/.Rprofile

ENTRYPOINT ["/bin/bash"]