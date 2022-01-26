```
docker run \
    -it \
    --rm \
    --user root \
    -p 8088:3838 \
    -v ${PWD}:/work \
    --name hello-testis \
    rnakato/singlecell_jupyter \
    R -e 'shiny::runApp("/work/hello-testis", host="0.0.0.0", port=3838)'
```