# scatac-mouse-testis

```
# processing
docker run \
    -it \
    --rm \
    --user root \
    -v ${PWD}:/work \
    --name scatac-mouse-testis \
    rnakato/singlecell_jupyter \
    /bin/bash
cd /work
Rscript R/010_qc.R
Rscript R/020_filtering.R
Rscript R/030_Integration_Un.R
Rscript R/040_Integration_Co.R
Rscript R/050_peakcall.R
```