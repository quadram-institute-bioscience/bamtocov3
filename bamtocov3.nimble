# Package

version       = "3.0.0"
author        = "Andrea Telatin"
description   = "BAM Coverage Toolkit"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.20", "docopt >= 0.6.8", "nim >= 1.6.6", "lapper >= 0.1.8", "taskpools >= 0.0.3", "colorize"

srcDir = "src"
binDir = "bin"

# src/main.nim is compiled to bin/bamto
namedBin["main"] = "bamto"

skipDirs = @["tests", "docs", "cli", "utils"]
skipFiles = @["example.bam"]

task test, "Build bamto and run Bats integration tests":
  exec "nimble build"
  exec "test/run.sh"
