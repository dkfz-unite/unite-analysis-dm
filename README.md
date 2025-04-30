# Differential Methylation Analysis (DM) Service

## General
DM analysis application wrapped with web API.

## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service:
- `UNITE_COMMAND` - command to run the analysis package (`Rscript`).
- `UNITE_COMMAND_ARGUMENS` - command arguments (`run.R {data}/{proc}_data.tsv {data}/{proc}_metadata.tsv {data}/{proc}_results.tsv`).
- `UNITE_SOURCE_PATH` - location of the source code in docker container (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/analysis`).
- `UNITE_LIMIT` - maximum number of concurrent jobs (`1` - process is heavy and uses a lot of CPU).

## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- Environment configuration and installation scripts: https://github.com/dkfz-unite/unite-env analysis service configuration and installation scripts: https://github.com/dkfz-unite/unite-env/tree/main/applications/unite-analysis-dm

### Docker
[Dockerfile](Dockerfile) is used to build an image of the application.
To build an image run the following command:
```
docker build -t unite.analysis.dm:latest .
```

All application components should run in the same docker network.
To create common docker network if not yet available run the following command:
```bash
docker network create unite
```


To run application in docker run the following command:
```bash
docker run \
--name unite.analysis.dm \
--restart unless-stopped \
--net unite \
--net-alias dm.analysis.unite.net \
-p 127.0.0.1:5302:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-v ./data:/mnt/data:rw \
-d \
unite.analysis.dm:latest
```

