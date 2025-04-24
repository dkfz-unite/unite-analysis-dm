FROM r-base:latest as base
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y r-base r-base-dev

FROM base as install
COPY ./src/install /src
WORKDIR /src

RUN apt-get update && apt-get install -y --allow-downgrades libcurl4  libcom-err2 libkrb5-dev libcurl4-openssl-dev libssl-dev libxml2-dev

RUN Rscript install.R
RUN apt-get clean

FROM install as final
COPY ./src/run /src
COPY ./app /app
WORKDIR /app
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
ENV ASPNETCORE_hostBuilder:reloadConfigOnChange=false
ENV UNITE_COMMAND="Rscript"
ENV UNITE_COMMAND_ARGUMENTS="run.R {data}/{proc}/metadata.tsv {data}/{proc}/options.json {data}/{proc}/results.csv.gz {data}/{proc}/results_reduced.csv.gz {data}/{proc}/results_annotated.csv.gz"
ENV UNITE_SOURCE_PATH="/src"
ENV UNITE_DATA_PATH="/mnt/data"
ENV UNITE_PROCESS_LIMIT="1"
EXPOSE 80
CMD ["/app/Unite.Commands.Web", "--urls", "http://0.0.0.0:80"]