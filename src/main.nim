## BamToCov 3 - entrypoint, compiled to bin/bamto.
##
## Dispatches to subcommands defined in src/cli/{name}.nim. Each subcommand
## exposes a `proc name(argv: var seq[string]): int`. Keep this file thin:
## the only job here is routing and the help/splash screen.

import std/[algorithm, sequtils, strutils, tables]
import colorize

import ./utils/utils

# Subcommands
import ./cli/cov
import ./cli/counts
import ./cli/contigs
import ./cli/convert
import ./cli/gen

# Dispatch table: subcommand name (and aliases) -> handler
type Handler = proc(argv: var seq[string]): int {.nimcall.}

let progs: Table[string, Handler] = {
  "cov":      Handler(cov.cov),
  "counts":   Handler(counts.counts),
  "cnt":      Handler(counts.counts),
  "contigs":  Handler(contigs.contigs),
  "convert":  Handler(convert.convert),
  "gen":      Handler(gen.gen),
}.toTable

# Short one-line description per primary subcommand, shown on the help screen
let helps = {
  "cov":      "per-base / per-target coverage (BED, BedGraph, WIG)",
  "counts":   "feature counts across one or more BAM files",
  "contigs":  "per-contig coverage statistics (CoverM-like)",
  "convert":  "convert annotation formats (GFF/GTF -> BED ...)",
  "gen":      "simulate / generate test data",
}.toTable

proc printHelp() =
  stderr.writeLine(splash().fgGreen)
  var names = toSeq(helps.keys)
  sort(names, system.cmp)
  for name in names:
    stderr.writeLine("  " & "· $1 $2".format(name.alignLeft(10), helps[name]))
  stderr.writeLine("""
Type 'bamto version' or 'bamto cite' to print the version and citation.
Add --help after each command to print its usage.""")

proc main(args: var seq[string]): int =
  if args.len == 1 and args[0] in ["version", "--version", "-v"]:
    echo version()
    return 0

  if args.len == 1 and args[0] in ["cite", "citation", "--cite", "--citation"]:
    echo "BamToCov ", version()
    echo "------------------------------------------------------------------------"
    echo "Birolo G, Telatin A."
    echo "BamToCov: an efficient toolkit for sequence coverage calculations."
    echo "Bioinformatics 2022, btac125. doi.org/10.1093/bioinformatics/btac125"
    echo ""
    echo "Repository:    https://github.com/telatin/bamtocov"
    return 0

  if args.len < 1 or args[0] notin progs:
    printHelp()
    if args.len > 0 and args[0] notin ["--help", "-h", "help"]:
      logErr("Unknown subcommand: " & args[0])
      return 1
    return 0

  # Valid subcommand: hand over the remaining arguments
  var pargs = args[1 .. ^1]
  return progs[args[0]](pargs)

when isMainModule:
  main_helper(main)
