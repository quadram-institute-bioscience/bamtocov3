## bamto gen - generate / simulate test data (reads, coverage, ...).

import docopt
import ../utils/utils

const usage = """
bamto gen: simulate test data

Usage:
  bamto gen [options]

Options:
  -h --help        Show this help
"""

proc gen*(argv: var seq[string]): int =
  let args = docopt(usage, argv = @["gen"] & argv, version = version())
  logErr("'bamto gen' is not implemented yet.")
  logDebug("args: ", $args)
  return 1
