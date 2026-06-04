## Coverage-specific helpers shared by the `cov` (and, later, `counts` /
## `contigs`) subcommands.
##
## This is the BamToCov 3 port of the original `covutils.nim`: genomic
## interval types, target (BED/GFF/GTF) parsing, and the sweep-line
## intersection machinery. Generic concerns (version, signal handling,
## logging) live in `utils.nim`, so they are intentionally absent here.

# Standard library
import std/[algorithm, strutils, tables]

# External dependencies
import hts

type
  region_t* = ref object
    chrom: string
    start: int
    stop: int
    name: string

type
  chrom_t* = int                  ## reference id (tid) inside the BAM header
  pos_t* = int64
  # Intervals carry a "label" with payload beyond their location. The plain
  # `interval_t` has no chromosome: it is used inside the target table, where
  # intervals are already grouped by chromosome.
  interval_t*[T]         = tuple[start, stop: pos_t, label: T]
  genomic_interval_t*[T] = tuple[chrom: chrom_t, start, stop: pos_t, label: T]
  target_index_t*        = tuple[chrom: chrom_t, interval: int]
  target_t*              = TableRef[chrom_t, seq[interval_t[string]]]
  raw_target_t*          = TableRef[string, seq[region_t]]

################################
# INTERVAL TYPES AND FUNCTIONS #
################################

proc is_null*(c: chrom_t): bool = c == -1

proc intersection_both*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): genomic_interval_t[tuple[l1: T1, l2: T2]] =
  (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), (i1.label, i2.label))

proc intersection_first*[T1, T2](i1: genomic_interval_t[T1], i2: genomic_interval_t[T2]): genomic_interval_t[T1] =
  if i1.chrom == i2.chrom:
    (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
  else:
    (i1.chrom, pos_t(0), pos_t(0), i1.label)

proc cookTarget*(orig: raw_target_t, bam: Bam): target_t =
  ## Resolve the chromosome names of a raw target against the BAM header,
  ## turning name-keyed regions into tid-keyed, sorted intervals.
  var chrom_map = newTable[string, chrom_t]()
  for t in bam.hdr.targets:
    chrom_map[t.name] = t.tid
  var cooked = newTable[chrom_t, seq[interval_t[string]]]()
  for chrom_str, intervals in orig:
    doAssert(not (":" in chrom_str), "bad target")
    if chrom_str notin chrom_map:
      raise newException(ValueError, "Target contig not found in BAM/CRAM header: " & chrom_str)
    let chrom = chrom_map[chrom_str]
    cooked[chrom] = @[]
    var last_start = 0
    for i in intervals:
      doAssert(i.chrom == chrom_str, "bad target")
      doAssert(i.start >= last_start)
      let name = if i.name == "": ($i.chrom & ":" & $i.start & "-" & $i.stop) else: i.name
      cooked[chrom].add((pos_t(i.start), pos_t(i.stop), name))
      last_start = i.start
  cooked

proc is_empty*[T](i: genomic_interval_t[T]): bool = i.start >= i.stop

# true when interval i1 is strictly before (no intersection) interval i2
proc `<<`*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i1.stop < i2.start
proc `<<`*[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i1.stop < i2.start
# true when interval i1 is strictly after (no intersection) interval i2
proc `>>`*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i2 << i1
proc `>>`*[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i2 << i1
# interval length
proc len*[T](i: genomic_interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
# true when interval is not empty
proc to_bool*[T](i: genomic_interval_t[T]): bool = i.start < i.stop

iterator intersections*[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): genomic_interval_t[tuple[l1: T, l2: string]] =
  ## Yield the intersection of `query` with each overlapping target interval.
  ## `idx` is carried between calls so the target is only ever scanned forward
  ## (the query stream is assumed coordinate-sorted).
  let chrom = query.chrom
  if chrom in target:
    let
      intervals = target[chrom]
      max_idx = len(intervals) - 1
    # if the chrom changes, start from the leftmost interval
    if idx.chrom != chrom:
      idx = (chrom, 0)
    # advance the target until it reaches the query interval
    while idx.interval < max_idx and intervals[idx.interval] << query:
      idx.interval += 1
    # yield all intersections
    var i = idx.interval
    while i <= max_idx and not (query << intervals[i]):
      yield intersection_both(query, intervals[i])
      i += 1

proc intersects*[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): bool =
  for i in intersections(query, target, idx):
    return to_bool(i)
  return false

# Accessors for the opaque region_t
proc start*(r: region_t): int {.inline.} = r.start
proc stop*(r: region_t): int {.inline.} = r.stop
proc chrom*(r: region_t): string {.inline.} = r.chrom
proc name*(r: region_t): string {.inline.} = r.name

#####################################
# TARGET (BED / GFF / GTF) PARSING  #
#####################################

# Converts a GTF line to a region object
proc gtf_line_to_region*(line: string, gffField = "exon", gffSeparator = ";", gffIdentifier = "gene_id"): region_t =
  var cse = line.strip().split('\t')

  if len(cse) < 8:
    stderr.write_line("[warning] skipping GTF line (fields not found):", line.strip())
    return nil

  # Skip non CDS fields (or the user-provided feature type)
  if cse[2] != gffField:
    return nil

  var
    s = parse_int(cse[3]) - 1
    e = parse_int(cse[4])
    reg = region_t(chrom: cse[0], start: s, stop: e)

  if len(cse) == 9:
    let gtfAnnotationFull = if gffSeparator in cse[8]: cse[8].split(gffSeparator)
                            else: @[cse[8]]
    for gffAnnotPartRaw in gtfAnnotationFull:
      let gffAnnotPart = gffAnnotPartRaw.strip(chars = {'"', '\'', ' '})
      if gffAnnotPart.startsWith(gffIdentifier):
        if "=" in gffAnnotPart:
          reg.name = gffAnnotPart.split("=")[1].strip(chars = {'"', '\'', ' '})
        elif " " in gffAnnotPart:
          reg.name = gffAnnotPart.split(" ")[1].strip(chars = {'"', '\'', ' '})
        else:
          reg.name = gffAnnotPart
        break
    if reg.name == "":
      reg.name = reg.chrom & ":" & $reg.start & "-" & $reg.stop
  return reg

# Converts a GFF line to a region object
proc gff_line_to_region*(line: string, gffField = "CDS", gffSeparator = ";", gffIdentifier = "ID"): region_t =
  var cse = line.strip().split('\t')

  # Skip unexpectedly short lines
  if len(cse) < 8:
    return nil

  # Skip non CDS fields (or the user-provided feature type)
  if cse[2] != gffField:
    return nil

  var s, e: int
  try:
    s = parse_int(cse[3]) - 1
    e = parse_int(cse[4])
  except:
    stderr.write_line("[warning] fields 4 and 5 are not integers):", cse[3], ", ", cse[4])
    return nil
  var reg = region_t(chrom: cse[0], start: s, stop: e)

  # In the future the 9th field could be required [TODO]
  if len(cse) == 9:
    try:
      for gffAnnotPartRaw in cse[8].split(gffSeparator):
        let gffAnnotPart = gffAnnotPartRaw.strip(chars = {'"', '\'', ' '})
        if gffAnnotPart.startsWith(gffIdentifier):
          let splittedField = gffAnnotPart.split("=")
          # Try splitting on "="
          if len(splittedField) == 2:
            reg.name = splittedField[1].strip(chars = {'"', '\'', ' '})
            break
          else:
            let resplittedField = gffAnnotPart.split(" ")
            if len(resplittedField) == 2:
              reg.name = resplittedField[1].strip(chars = {'"', '\'', ' '})
              break
            else:
              reg.name = "Error"
              break
    except Exception as e:
      stderr.write_line("[warning] fields 8 is not a string):", cse[8], "\n  ", e.msg)
      return nil

  return reg

# Converts a BED line to a region object.
#
# Note: GFF/GTF input is NOT auto-detected from the *content* of a line here.
# The format is decided upstream by file extension or `--target-format`
# (see `cov.nim`), and only then is `gff_to_table` / `gtf_to_table` used.
# A previous attempt to sniff GFF inside this proc (`if len(cse) == 9`) was
# dead code: `split('\t', 5)` performs at most 5 splits, so `cse` never holds
# more than 6 fields and the 9-field branch could never fire. If a GFF/GTF
# line nonetheless reaches this proc (e.g. the user forced `--target-format
# bed`), its 2nd/3rd columns are non-numeric and the coordinate parse below
# would raise an unhandled `ValueError`. We catch that and emit an actionable
# message instead of crashing.
proc bed_line_to_region*(line: string): region_t =
  var cse = line.strip().split('\t', 5)

  if len(cse) < 3:
    stderr.write_line("[warning] skipping bad bed line:", line.strip())
    return nil

  var s, e: int
  try:
    s = parse_int(cse[1])
    e = parse_int(cse[2])
  except ValueError:
    stderr.write_line("[warning] columns 2 and 3 are not integers; this does not look like BED. ",
                      "For GFF/GTF use a .gff/.gtf extension or pass --target-format. Skipping line: ",
                      line.strip())
    return nil

  var reg = region_t(chrom: cse[0], start: s, stop: e)
  if len(cse) > 3:
    reg.name = cse[3]
  return reg

# Read a delimited target file (BED/GFF/GTF) into name-keyed, per-chrom regions.
proc bed_to_table*(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  if bed == "nil":
    return bed_regions

  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, sort each contig's intervals by start
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc gtf_to_table*(bed: string, gffField, gffSeparator, gffIdentifier: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("##FASTA"):
      break
    if $kstr.s[0] == "#":
      continue
    var v = gtf_line_to_region($kstr.s, gffField, gffSeparator, gffIdentifier)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc gff_to_table*(bed: string, gffField, gffSeparator, gffIdentifier: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("##FASTA"):
      break
    if $kstr.s[0] == "#":
      continue

    var v: region_t
    try:
      v = gff_line_to_region($kstr.s, gffField, gffSeparator, gffIdentifier)
    except Exception as e:
      stderr.write_line("[GFF/GTF error]:", e.msg)
      continue

    if v == nil:
      continue

    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc target_names_in_order*(path: string, formatGff, formatGtf: bool, gffField, gffSeparator, gffIdentifier: string): seq[string] =
  ## Return the target interval names in the order they appear in the file,
  ## used to lay out the report table rows deterministically.
  if path == "nil":
    return @[]

  var hf = hts.hts_open(cstring(path), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil

  while hts_getline(hf, cint(10), addr kstr) > 0:
    let line = $kstr.s
    if formatGff or formatGtf:
      if line.startsWith("##FASTA"):
        break
      if line.len == 0 or line[0] == '#':
        continue
    else:
      if line.startsWith("track "):
        continue
      if line.len == 0 or line[0] == '#':
        continue

    let region =
      if formatGff:
        gff_line_to_region(line, gffField, gffSeparator, gffIdentifier)
      elif formatGtf:
        gtf_line_to_region(line, gffField, gffSeparator, gffIdentifier)
      else:
        bed_line_to_region(line)

    if region != nil:
      result.add(if region.name == "": region.chrom & ":" & $region.start & "-" & $region.stop else: region.name)

  hts.free(kstr.s)
