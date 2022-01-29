markerGenes  <- c(
    "Foxo1", "Neurog3", "Tspan8", "Sall4", "Sohlh1", "Sohlh2",
    "Gfra1", "Id4", "Nanos2", "Zbtb16",  # undifferentiated SG
    "Etv5", "Dmrt1", "Kit", "Rhox13", "Stra8", "Uchl1", # SPG
    "Prss50", "Tex101", "Ly6k",  # preLeptotene
    "Spag6", "Sycp3", "Dmc1", "Piwil1",  "Id4", "Pgk2",  "Spag6", "Tbpl1",  # Spermatocytes
    "Pou5f2", "Kctd9", "Mybl1", "Adam3",  # Diplotene
    "Acrv1", "Spaca1", "Tsga8", "Tex21",  # Spermatids
    "Gapdhs", "Tnp1", "Tnp2", "Prm1", "Prm2",  # Elongating
    "Clu", "Ctsl", "Sox9", "Amh", "Cldn11",  # Sertoli
    "Cyp17a1", "Cyp11a1", "Star", "Fabp3",  # Leydig
    "Pdgfra",  # SLC
    "Nr5a1", "Acta2", "Igf1"
  )

write.table(
    markerGenes, 
    file = "gene_sets/markerGenes.txt", 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE)

markerGenes <- read.table("gene_sets/markerGenes.txt")[["V1"]]