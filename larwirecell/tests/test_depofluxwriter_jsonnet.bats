@test "compile wct configuration" {

    local name="test_depofluxwriter"
    local mydir="$(dirname "$BATS_TEST_FILENAME")"
    local cfg=$mydir/depofluxwriter/cfg

    # these are bogus parameters
    run jsonnet \
        --ext-code G4RefTime=0.0 \
        --ext-code lifetime=1.0 \
        --ext-code DT=1.0 \
        --ext-code DL=1.0 \
        --ext-code driftSpeed=1.0 \
        -o $name.json \
        $cfg/$name.jsonnet
    echo "$output"
    [[ "$status" -eq 0 ]]
}
