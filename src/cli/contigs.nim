## bamto contigs - per-contig coverage statistics (CoverM-like).

import docopt
import ../utils/utils

const usage = """
bamto contigs: per-contig coverage statistics

Usage:
  bamto contigs [options] <input.bam>...

Options:
  -h --help        Show this help
"""

proc contigs*(argv: var seq[string]): int =
  let args = docopt(usage, argv = @["contigs"] & argv, version = version())
  logErr("'bamto contigs' is not implemented yet.")
  logDebug("args: ", $args)
  return 1
