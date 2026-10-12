# Differential Methylation Analysis (DM) Service

Shared architecture, development workflows, and agent instructions:
[UNITE project documentation](https://github.com/dkfz-unite/.github/blob/main/README.md).

## General
Differential methylation analysis of Illumina methylation array data with [minfi](https://bioconductor.org/packages/minfi) and [limma](https://bioconductor.org/packages/limma), wrapped with web API.

The R script ([src/run/run.R](src/run/run.R)) is run by the [commands](https://github.com/dkfz-unite/unite-commands) web service, which runs a configured command for each request. The service is used by the [analysis service](https://github.com/dkfz-unite/unite-analysis) for [DM analysis](https://github.com/dkfz-unite/unite-analysis/blob/main/Docs/api-analysis-dm.md): the analysis service downloads the IDAT files and writes the input files to the shared analysis data directory, and calls the service with the analysis id as process key.

The script performs the following steps:
1. Reads the IDAT files of the samples (`read.metharray`) and preprocesses them with the selected minfi method.
2. Calculates M-values (`getM`).
3. Compares the two conditions with `limma` (`lmFit`, `eBayes`, `topTable` for the condition coefficient). The condition of the last sample in the metadata is the reference.
4. Rounds `logFC` and `adj.P.Val` to 4 decimals, removes CpGs with non-finite `logFC`, and keeps a random sample of at most 50 000 CpGs in each of 10 `logFC` bins (this is done twice: before and after annotation).
5. Annotates the CpGs with the Illumina EPIC annotation (`IlluminaHumanMethylationEPICanno.ilm10b4.hg19`).
6. Writes the results.

The script uses the `IlluminaHumanMethylationEPICmanifest` manifest. R packages: `minfi`, `limma`, `jsonlite`, `dplyr` and the EPIC manifest and annotation packages. [src/install/install.R](src/install/install.R) installs only `dplyr`; the other packages are expected from the base image `ghcr.io/dkfz-unite/docker-rdev-minfi` (their install lines are commented out).

## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service. The [Dockerfile](Dockerfile) sets the following values:
- `UNITE_COMMAND` - command to run the analysis script (`Rscript`).
- `UNITE_COMMAND_ARGUMENTS` - command arguments (`run.R {data}/{proc}/metadata.tsv {data}/{proc}/options.json {data}/{proc}/results.tsv`).
  - `{data}` is replaced with `UNITE_DATA_PATH`, `{proc}` with the process key (analysis id).
- `UNITE_SOURCE_PATH` - location of the script in docker container, the working directory of the command (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/data`).
- `UNITE_PROCESS_LIMIT` - maximum number of concurrent processes (`1` - process is heavy and uses a lot of CPU and memory). If the limit is reached, the service returns `501`.

## Usage
- Place `metadata.tsv`, `options.json` and the IDAT files to the `{proc}` subdirectory of the data directory (`/mnt/data/{proc}` in the container).
- Send a `POST` request to `/api/run?key={proc}`.
- The service runs `Rscript run.R /mnt/data/{proc}/metadata.tsv /mnt/data/{proc}/options.json /mnt/data/{proc}/results.tsv` and returns `200` with the script output, or `500` if the script fails.

### Input files
`metadata.tsv` - samples:
- `sample_id` - sample identifier.
- `path` - path of the sample IDAT files without the channel suffix and extension (`{path}_Red.idat` and `{path}_Grn.idat` are read). The path has to be valid in the runner container.
- `condition` - condition of the sample. There should be exactly two conditions.

```tsv
sample_id	path	condition
D001	/mnt/data/{proc}/D001/203711230017_R01C01	8b4e6a2c-2b4b-4e3e-9c6d-5a4f1d3b9e0f
D002	/mnt/data/{proc}/D002/203711230017_R02C01	2e6ee2ff-d96d-422d-8a0d-a49a412baa5c
```

The analysis service stores the IDAT files of each sample in the `{proc}/{sample_id}` directory, writes absolute paths (so the data directory has to be mounted at the same path as in the analysis service, `/mnt/data`) and uses the dataset id as condition.

`options.json` - preprocessing options:
- `pp` - minfi preprocessing method: `"Illumina"` (`preprocessIllumina`, also used for any other value), `"SWAN"` (`preprocessSWAN`), `"Quantile"` (`preprocessQuantile`), `"Noob"` (`preprocessNoob`), `"Raw"` (`preprocessRaw`).

```json
{
  "pp": "Illumina"
}
```

### Output files
`results.tsv` - annotated results (subsampled, see above):
- `CpgId` - CpG probe id.
- `logFC` - difference of mean M-values (test condition minus reference condition).
- `adj.P.Val` - adjusted p-value.
- `UCSC_RefGene_Name`, `Regulatory_Feature_Name` - annotation columns.
- Annotation columns with `Enhancer` in their name (e.g. `Phantom4_Enhancers`, `Phantom5_Enhancers`).
- `logFc` - copy of `logFC`.
- `logFc_bin` - `logFC` bin used for subsampling.

The analysis service reads `CpgId`, `logFc`, `adj.P.Val`, `Regulatory_Feature_Name`, `Phantom4_Enhancers`, `Phantom5_Enhancers` and `UCSC_RefGene_Name` to create a binned volcano plot file.

## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- [Environment configuration and installation scripts](https://github.com/dkfz-unite/unite-env)
- [DM analysis service configuration and installation scripts](https://github.com/dkfz-unite/unite-env/tree/main/applications/unite-analysis-dm)

### Docker
The image of the service is available in our [registry](https://github.com/dkfz-unite/unite-analysis-dm/pkgs/container/unite-analysis-dm).

The image is built by the [docker-images](.github/workflows/docker-images.yml) workflow on pushes to `main`: it checks out [unite-commands](https://github.com/dkfz-unite/unite-commands), publishes it with `publish-linux-x64.sh` and `publish-linux-arm64.sh`, copies the result (`commands` executable) to `./app`, builds amd64 and arm64 images with the [Dockerfile](Dockerfile) and pushes them under one `latest` manifest.

To build an image locally, publish unite-commands for the target platform to `./app` first (the `app` directory in the repository is replaced by the workflow and may not match the target platform or the current unite-commands version):
```bash
# in the unite-commands repository
sh publish-linux-x64.sh # or publish-linux-arm64.sh

# in this repository
rm -rf app
cp -r ../unite-commands/publish/linux-x64 app
chmod a+rx app/commands
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
-p 127.0.0.1:5306:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-v unite.analysis.data:/mnt/data:rw \
-d \
unite.analysis.dm:latest
```

The volume `unite.analysis.data` is the analysis data directory shared with the analysis service and has to be mounted at the same path (`/mnt/data`), because `metadata.tsv` contains absolute paths.
