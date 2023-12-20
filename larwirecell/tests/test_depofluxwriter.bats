#!/usr/bin/env bats

function cd_tmp () {
    if [[ -n "$WCT_BATS_TMPDIR" ]] ; then
        mkdir -p "$WCT_BATS_TMPDIR"
        cd "$WCT_BATS_TMPDIR"
        return
    fi
    cd "$BATS_TEST_TMPDIR"
}

setup_file () {
    local mydir="$(dirname "$BATS_TEST_FILENAME")"
    local name="$(basename "$BATS_TEST_FILENAME" .bats)"

    export FHICL_FILE_PATH="$mydir/depofluxwriter/fcl:$FHICL_FILE_PATH"
    export WIRECELL_PATH="$mydir/depofluxwriter/cfg:$WIRECELL_PATH"

    cd_tmp

    fcl="${name}.fcl"
    out="${name}_artroot.root"
    if [[ -f "$out" ]] ; then
        echo "Existing artroot file: $out" 1>&3
    else
        echo art -n 1 -o "$out" -c "$fcl" 1>&3
        art -n 1 -o "$out" -c "$fcl"
        echo "$output"
        [[ "$status" -eq 0 ]]
    fi

    local data="${name}_artroot.root"
    local hist="${name}_hist.root"
    local plot="${name}_hist.pdf"

    if [[ -f "$hist" ]] ; then
        echo "Already have $hist" 1>&3
    else
        root -b -q $mydir/depofluxwriter/root/collect.C'("'"$data"'","'"$hist"'")'
    fi
    if [[ -f "$plot" ]] ; then
        echo "Already have $plot" 1>&3
    else
        root -b -q $mydir/depofluxwriter/root/plot.C'("'"$hist"'","'"$plot"'")'
    fi
}

@test "Assure SimEnergyDeposits match" {
    cd_tmp

    local src=wcls-sim-drift-deposource.log
    local tgt=wcls-sim-drift-depowriter.log
    run diff -u \
     <(grep '^SimDepoSetSource IonAndScint sed' $src | sed 's,SimDepoSetSource IonAndScint ,,' ) \
     <(grep '^DepoFluxWriter IonAndScint sed'   $tgt | sed 's,DepoFluxWriter IonAndScint ,,' ) 
    echo "$output"
    [[ "$status" -eq 0 ]]

}

# Depo columns
#  1 prefix
#  2 label
#  3 "depo"
#  4 child=0/parent=1
#  5 ndepos
#  6 index in array
#  7 parent ID
#  8 ID
#  9 PDG
# 10 charge
# 11 energy
# 12 time
# 13 x
# 14 y
# 15 z
# 16 longitudinal extent
# 17 transverse extent


@test "Assure IDepos match" {
    cd_tmp

    local src=wcls-sim-drift-deposource.log
    local tgt=wcls-sim-drift-depowriter.log

    run diff  \
     <(grep '^SimDepoSetSource IonAndScint depo 0 ' $src| awk '{print $8}') \
     <(grep '^DepoFluxWriter simpleSC depo 1 '     $tgt | awk '{print $8}')
    echo "$output"

    # Do not allow change in order or additions.  Omissions are okay.
    [[ -z "$(echo $output | grep '>')" ]]

}
