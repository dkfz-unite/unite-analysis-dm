FROM ubuntu:latest AS base
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y r-base r-base-dev

FROM base AS install
COPY ./src/install /src
WORKDIR /src

RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev

RUN Rscript install.R
RUN apt-get clean

FROM install AS final
COPY ./src/run /src
COPY ./app /app
WORKDIR /app
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
ENV ASPNETCORE_hostBuilder:reloadConfigOnChange=false
ENV UNITE_COMMAND="Rscript"
ENV UNITE_COMMAND_ARGUMENTS="run.R {data}/{proc}/metadata.tsv {data}/{proc}/options.json {data}/{proc}/results.tsv.gz"
ENV UNITE_SOURCE_PATH="/src"
ENV UNITE_DATA_PATH="/mnt/data"
ENV UNITE_PROCESS_LIMIT="1"
EXPOSE 80
CMD ["/app/commands", "--urls", "http://0.0.0.0:80"]