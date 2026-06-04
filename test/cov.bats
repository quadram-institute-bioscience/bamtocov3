#!/usr/bin/env bats

setup() {
    load 'test_helper/common-setup'
    _common_setup
}

@test "cov emits BED coverage for sorted BAM" {
    run "$BIN" cov "$DATA_DIR/cov/mini-sorted.bam"

    assert_success
    assert_output --partial $'EMPTY\t0\t414\t0'
    assert_output --partial $'NC_001416.1\t99\t199\t1'
    assert_output --partial $'NC_001416.1\t1299\t2800\t0'
}

@test "cov writes coverage to an output file" {
    coverage="$BATS_TEST_TMPDIR/mini.bedgraph"

    run "$BIN" cov --output "$coverage" "$DATA_DIR/cov/mini-sorted.bam"
    assert_success

    run test -s "$coverage"
    assert_success

    run grep -q $'NC_001416.1\t99\t199\t1' "$coverage"
    assert_success
}

@test "cov writes a summary when coverage output is skipped" {
    summary="$BATS_TEST_TMPDIR/mini-summary.tsv"

    run "$BIN" cov --summary-only --summary "$summary" "$DATA_DIR/cov/mini-sorted.bam"
    assert_success

    run test -s "$summary"
    assert_success

    run grep -q "^interval" "$summary"
    assert_success

    run grep -q "^NC_001416.1" "$summary"
    assert_success
}

@test "cov emits fixed-step WIG with explicit format and span" {
    run "$BIN" cov --format wig --span 100 "$DATA_DIR/cov/mini-sorted.bam"

    assert_success
    assert_output --partial "fixedStep chrom=EMPTY start=1 step=100 span=100"
    assert_output --partial "fixedStep chrom=NC_001416.1 start=1 step=100 span=100"
}

@test "cov accepts the camelCase extendReads alias" {
    run "$BIN" cov --extendReads 100 "$DATA_DIR/cov/mini-sorted.bam"

    assert_success
}

@test "cov fails cleanly for a missing BAM" {
    run "$BIN" cov "$DATA_DIR/cov/does-not-exist.bam"

    assert_failure 1
    assert_output --partial "Input BAM file"
    assert_output --partial "not found"
}
