#!/usr/bin/env bats

setup() {
    load 'test_helper/common-setup'
    _common_setup
}

@test "bamto binary is available" {
    [ -x "$BIN" ]
}

@test "bamto binary matches the current shell" {
    run _assert_bamto_binary_matches_current_shell

    assert_success
}

@test "version command reports package version" {
    run "$BIN" version

    assert_success
    assert_output "3.0.0"
}

@test "top-level help advertises the subcommands" {
    run "$BIN" --help

    assert_success
    assert_output --partial "BamToCov 3.0.0"
    assert_output --partial "per-base / per-target coverage"
    assert_output --partial "feature counts across one or more BAM files"
    assert_output --partial "per-contig coverage statistics"
    assert_output --partial "convert annotation formats"
    assert_output --partial "simulate / generate test data"
}

@test "cov help uses the new subcommand spelling" {
    run "$BIN" cov --help

    assert_success
    assert_output --partial "bamto cov [options] [<bam>...]"
    refute_output --partial "Usage: bamtocov"
}

@test "cov help advertises the reshaped v3 options" {
    run "$BIN" cov --help

    assert_success
    assert_output --partial "--min-mapq <N>"
    assert_output --partial "--exclude-flag <FLAG>"
    assert_output --partial "--extend-reads <BP>"
    assert_output --partial "--target-format <FMT>"
    assert_output --partial "--gff-feature <TYPE>"
    assert_output --partial "--gff-attribute <KEY>"
    assert_output --partial "--format <FMT>"
    assert_output --partial "--output <FILE>"
    assert_output --partial "--summary <FILE>"
    assert_output --partial "--summary-only"
    assert_output --partial "--low-cov <N>"
    assert_output --partial "--span <BP>"
    assert_output --partial "--span-op <FUNC>"
}

@test "cov help no longer advertises old bamtocov option names" {
    run "$BIN" cov --help

    assert_success
    refute_output --partial "--report "
    refute_output --partial "--skip-output"
    refute_output --partial "--wig"
    refute_output --partial "--op "
    refute_output --partial "--mapq"
    refute_output --partial "--extendReads"
}

@test "removed cov options fail with guidance" {
    run "$BIN" cov --report "$BATS_TEST_TMPDIR/old.tsv" "$DATA_DIR/cov/mini-sorted.bam"

    assert_failure 1
    assert_output --partial "--report was removed"
    assert_output --partial "--summary"
}

@test "counts alias cnt routes to counts" {
    run "$BIN" cnt --help

    assert_success
    assert_output --partial "Usage:"
    assert_output --partial "bamto counts"
}

@test "unknown subcommand fails" {
    run "$BIN" nope

    assert_failure 1
    assert_output --partial "Unknown subcommand: nope"
}
