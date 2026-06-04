## bamto counts - feature counts across multiple BAM files.
## (was: the original `bamtocounts`)

import docopt
import ../utils/utils

const usage = """
bamto counts: feature counts from one or more BAM files

Usage:
  bamto counts [options] <annotation> <input.bam>...

Options:
  -h --help        Show this help
"""

proc counts*(argv: var seq[string]): int =
  let args = docopt(usage, argv = @["counts"] & argv, version = version())
  logErr("'bamto counts' is not implemented yet.")
  logDebug("args: ", $args)
  return 1
