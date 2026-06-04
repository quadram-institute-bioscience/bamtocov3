#!/usr/bin/env bash

_common_setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'

    PROJECT_ROOT="$(cd "$(dirname "$BATS_TEST_FILENAME")/.." >/dev/null 2>&1 && pwd)"
    BIN="$PROJECT_ROOT/bin/bamto"
    DATA_DIR="$PROJECT_ROOT/test/data"

    export PROJECT_ROOT
    export BIN
    export DATA_DIR
}

_assert_bamto_binary_matches_current_shell() {
    local host_os host_arch bin_info run_error

    host_os="$(uname -s)"
    host_arch="$(uname -m)"

    if command -v file >/dev/null 2>&1; then
        bin_info="$(file "$BIN")"
        case "$host_os" in
            Darwin)
                if [[ "$bin_info" != *"Mach-O"* ]]; then
                    echo "Expected a macOS Mach-O binary for $host_os/$host_arch"
                    echo "Got: $bin_info"
                    return 1
                fi
                ;;
            Linux)
                if [[ "$bin_info" != *"ELF"* ]]; then
                    echo "Expected a Linux ELF binary for $host_os/$host_arch"
                    echo "Got: $bin_info"
                    return 1
                fi
                ;;
            *)
                echo "Unsupported test host: $host_os/$host_arch"
                echo "Binary: $bin_info"
                return 1
                ;;
        esac
    else
        bin_info="file(1) not available"
    fi

    if ! run_error="$("$BIN" version 2>&1 >/dev/null)"; then
        echo "Current shell cannot execute $BIN"
        echo "Host: $host_os/$host_arch"
        echo "Binary: $bin_info"
        if [[ -n "$run_error" ]]; then
            echo "$run_error"
        fi
        return 1
    fi
}
