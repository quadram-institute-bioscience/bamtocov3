## bamto convert - convert between annotation formats (GFF/GTF -> BED, ...).

import docopt
import ../utils/utils

const usage = """
bamto convert: convert between annotation formats

Usage:
  bamto convert [options] <input>

Options:
  -h --help        Show this help
"""

proc convert*(argv: var seq[string]): int =
  let args = docopt(usage, argv = @["convert"] & argv, version = version())
  logErr("'bamto convert' is not implemented yet.")
  logDebug("args: ", $args)
  return 1
