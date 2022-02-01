# sh ../scripts/annotateEnhancers_distal.sh > ../log/annotateEnhancers_distal.log 2>&1
PARDIR="Results/CREs/exp_07"
OUT_META="Results/CREs/exp_07/metaplot"

mkdir -p $OUT_META/{scATAC,H3K27ac}

for i in {1..10}
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
