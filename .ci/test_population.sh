export BASE=$DATADIR/cellular/circuit-1k

srun appositionizer \
    --output=$PWD \
    --from head \
    --to tail \
    --appositionspace=axodendritic \
    --recipe $BASE/bioname/recipe.json \
    --circuit-config $BASE/circuit_config_multipop.json |& tee output

set +x

# Output for 550 edeges between the (original 1k) ids <50 and >=900
grep "Number of edeges written: 550" output &> /dev/null || exit 1
