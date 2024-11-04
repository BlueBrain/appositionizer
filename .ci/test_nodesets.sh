export BASE=$DATADIR/cellular/circuit-1k

srun appositionizer \
    --output=$PWD \
    --from-nodeset test \
    --to-nodeset test \
    --recipe $BASE/bioname/recipe.json \
    --circuit-config $BASE/circuit_config.json |& tee output

set +x

grep "Number of cells written: 50" output &> /dev/null || exit 1
