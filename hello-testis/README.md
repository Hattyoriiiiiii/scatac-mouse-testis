```
docker run \
    -it \
    --user root \
    -p 8080:3838 \
    -v ${PWD}:/work \
    --name hello-testis \
    rnakato/singlecell_jupyter \
    Rscript -e /work/hello-testis/app.R
```