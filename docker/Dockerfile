FROM debian:stable
RUN apt update
RUN apt install -y g++ libgdal-dev libshp-dev libboost-all-dev gdal-bin cmake
WORKDIR PHES_Searching
COPY . .
RUN ./build.sh
