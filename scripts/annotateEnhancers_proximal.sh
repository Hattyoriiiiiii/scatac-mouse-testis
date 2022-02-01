# sh ../scripts/annotateEnhancers_proximal.sh > ../log/annotateEnhancers_proximal.log 2>&1
PARDIR="results/07_p2g_proximal/CREs/exp_01"
OUT_META="results/07_p2g_proximal/CREs/exp_01/metaplot"

mkdir -p $OUT_META/{scATAC,H3K27ac}

for i in {1..11}
do
    cluster_i=$(printf "%02d\n" ${i})
    echo $cluster_i

    ### HOMER ---------------
    findMotifsGenome.pl \
        ${PARDIR}/bed/cluster_${cluster_i}.bed \
        mm10 \
        ${PARDIR}/motif/cluster_${cluster_i}/ \
        -size given &

    ### deepTools ---------------
    ##### H3K27ac
    computeMatrix reference-point \
        --referencePoint center \
        -a 3000 -b 3000 \
        -R ${PARDIR}/bed/cluster_${cluster_i}.bed \
        -S Inputs/bigwig/H3K27ac/*.normalized.bw \
        --skipZeros \
        -p 6 \
        -o ${OUT_META}/H3K27ac/H3K27ac_level_c${cluster_i}.norm.gz && \
    plotProfile \
        -m ${OUT_META}/H3K27ac/H3K27ac_level_c${cluster_i}.norm.gz \
        -out ${OUT_META}/H3K27ac/H3K27ac_profile_c${cluster_i}.norm.pdf \
        --plotHeight 14 \
        --plotWidth 22 \
        --perGroup  --plotTitle "" \
        -T "P2G sites"  -z "" \
        --startLabel "" \
        --endLabel "" &

    ##### scATAC
    computeMatrix reference-point \
        --referencePoint center \
        -a 3000 -b 3000 \
        -R ${PARDIR}/bed/cluster_${cluster_i}.bed \
        -S Inputs/bigwig/scATAC/*.bw \
        --skipZeros \
        -p 6 \
        -o ${OUT_META}/scATAC/scATAC_level_c${cluster_i}.norm.gz && \
    plotProfile \
        -m ${OUT_META}/scATAC/scATAC_level_c${cluster_i}.norm.gz \
        -out ${OUT_META}/scATAC/scATAC_profile_c${cluster_i}.norm.pdf \
        --plotHeight 14 \
        --plotWidth 22 \
        --perGroup  --plotTitle "" \
        -T "P2G sites"  -z "" \
        --startLabel "" \
        --endLabel "" &
done
