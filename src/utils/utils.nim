## Generic utilities shared across all BamToCov3 subcommands.
##
## Program-specific helpers live in src/utils/{program_name}.nim;
## anything generic (version, logging, the splash screen, the main
## dispatch helper) belongs here.

import std/[os, strutils, sugar]
when not defined(windows):
  import std/posix

# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------

const NimblePkgVersion {.strdefine.} = "<NimblePkgVersion>"

proc version*(): string =
  ## Package version, injected at compile time by nimble.
  if NimblePkgVersion.len == 0 or NimblePkgVersion == "<NimblePkgVersion>":
    "0.0.0"
  else:
    NimblePkgVersion

# ---------------------------------------------------------------------------
# Global switches (driven by environment variables)
# ---------------------------------------------------------------------------

let
  debug*: bool = getEnv("BAMTO_DEBUG", "0") notin ["", "0"]
  quiet*: bool = getEnv("BAMTO_QUIET", "0") notin ["", "0"]

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

proc logMsg*(things: varargs[string, `$`]) =
  ## Informational message to stderr (suppressed when BAMTO_QUIET is set).
  if quiet: return
  stderr.writeLine(things.join(" "))

proc logErr*(things: varargs[string, `$`]) =
  ## Error message to stderr (always shown).
  stderr.writeLine("[error] " & things.join(" "))

proc logDebug*(things: varargs[string, `$`]) =
  ## Debug message to stderr (only when BAMTO_DEBUG is set).
  if debug:
    stderr.writeLine("[debug] " & things.join(" "))

# ---------------------------------------------------------------------------
# Splash screen
# ---------------------------------------------------------------------------

const Logo* = """
  ___           _____     ___           ____
 | _ ) __ _ _ _|_   _|__ / __|_____ __ |__ /
 | _ \/ _` | '  \| |/ _ \ (__/ _ \ V /  |_ \
 |___/\__,_|_|_|_|_|\___/\___\___/\_/  |___/
"""

proc splash*(): string =
  ## The ASCII logo plus the version line.
  Logo & "\n BamToCov " & version() & " - BAM Coverage Toolkit\n"

# ---------------------------------------------------------------------------
# main_helper: wraps the dispatcher with signal/exception handling
# ---------------------------------------------------------------------------

proc main_helper*(main_func: var seq[string] -> int) =
  ## Run `main_func` with the command-line arguments, handling broken pipes,
  ## Ctrl-C and uncaught exceptions cleanly.
  var args: seq[string] = commandLineParams()
  when defined(windows):
    try:
      quit(main_func(args))
    except IOError:
      quit(0)
    except Exception:
      stderr.writeLine(getCurrentExceptionMsg())
      quit(2)
  else:
    signal(SIGPIPE, cast[typeof(SIG_IGN)](proc(sig: cint) =
      logDebug("handled SIGPIPE")
      quit(0)
    ))

    proc handler() {.noconv.} =
      if not quiet:
        stderr.writeLine("[Quitting on Ctrl-C]")
      quit(1)
    setControlCHook(handler)

    try:
      let exitStatus = main_func(args)
      logDebug("exiting ", $exitStatus)
      quit(exitStatus)
    except IOError:
      logDebug("IOError (broken pipe)")
      quit(1)
    except Exception:
      stderr.writeLine(getCurrentExceptionMsg())
      quit(2)
