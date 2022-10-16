```
# at scatac-mouse-testis
docker run \
    -it \
    --rm \
    --user root \
    -p 58080:3838 \
    -v ${PWD}:/work \
    --name hello-testis \
    rnakato/singlecell_jupyter \
    R -e 'shiny::runApp("/work/hello-testis", host="0.0.0.0", port=3838)'

# temp run
docker run \
    -it \
    --rm \
    --user root \
    -p 80:3838 \
    -v ${PWD}:/work \
    --name hello-testis \
    hattyoriiiiiii/single-cell-analysis:1.0.0 \
    /bin/bash

R
BiocManager::install("ChIPpeakAnno")
install.packages("shinybusy")
shiny::runApp("/work/hello-testis", host="0.0.0.0", port=3838)
```

```
docker run \
    -it \
    --rm \
    --user root \
    -p 80:3838 \
    -v ${PWD}:/work \
    --name hello-testis \
    hattyoriiiiiii/single-cell-analysis:1.0.0 \
    R -e 'shiny::runApp("/work/hello-testis", host="0.0.0.0", port=3838)'
```