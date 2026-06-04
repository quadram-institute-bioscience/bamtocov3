## bamto cov - per-base / per-target coverage from a BAM/CRAM file.
##
## This is the BamToCov 3 home of the original `bamtocov` tool: a sweep-line
## coverage calculator that emits BED / BedGraph / WIG, optionally restricted
## to a target (BED/GFF/GTF) and/or summarized into a per-interval report.

# Standard library
import std/[algorithm, heapqueue, os, sequtils, sets, strutils, tables]

# External dependencies
import docopt, hts

# Local modules
import ../utils/utils
import ../utils/coverage_utils

type
  coverage_t = object
    forward: int
    reverse: int

var debugCov = false

proc db(things: varargs[string, `$`]) =
  stderr.write("[debug]")
  for t in things:
    if t.len > 0 and t[0] == ',':
      stderr.write(t)
    else:
      stderr.write(" " & t)
  stderr.write("\n")

template dbEcho(things: varargs[string, `$`]) =
  if debugCov:
    db(things)

template covAssert(condition: bool, message: string) =
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1)

proc `$`[T](i: genomic_interval_t[T]): string =
  $i.chrom & ":" & $i.start & "-" & $i.stop & $i.label

type
  input_option_t = tuple[min_mapping_quality: uint8, eflag: uint16, physical: bool, extendFrag: int, target: raw_target_t]

proc alignment_stream(bam: Bam, opts: input_option_t, target: target_t): iterator (): genomic_interval_t[bool] =
  result = iterator(): genomic_interval_t[bool] {.closure.} =
    var
      o = opts
      target_idx: target_index_t
    for r in bam:
      # alignment filter
      if r.mapping_quality < o.min_mapping_quality or (r.flag and o.eflag) != 0:
        continue

      # alignment processing
      var stop: pos_t = 0
      if o.physical:
        if r.isize > 0:
          stop = r.start + r.isize
        else:
          continue # skip the mate with negative insert size
      else:
        stop = r.stop
      let i = (r.tid, pos_t(r.start), stop, r.flag.reverse)

      # return alignment if it intersects the target (or if there is no target)
      if len(target) == 0 or i.intersects(target, target_idx):
        yield i

# COVERAGE FUNCTIONS #
proc newCov(f = 0, r = 0): coverage_t =
  coverage_t(forward: f, reverse: r)

proc inc(c: var coverage_t, reverse = false) =
  if reverse == false:
    c.forward += 1
  else:
    c.reverse += 1

proc dec(c: var coverage_t, reverse = false) =
  if reverse == false:
    c.forward -= 1
  else:
    c.reverse -= 1

proc tot(c: coverage_t): int =
  c.forward + c.reverse

proc max(c1: coverage_t, c2: coverage_t): coverage_t =
  newCov(max(c1.forward, c2.forward), max(c1.reverse, c2.reverse))
proc min(c1: coverage_t, c2: coverage_t): coverage_t =
  newCov(min(c1.forward, c2.forward), min(c1.reverse, c2.reverse))
proc `+`(c1: coverage_t, c2: coverage_t): coverage_t =
  newCov(c1.forward + c2.forward, c1.reverse + c2.reverse)

proc `/`(c: coverage_t, by: float): tuple[forward: float, reverse: float] =
  (float(c.forward)/by, float(c.reverse)/by)
proc `*`(c: coverage_t, by: int): coverage_t =
  newCov(c.forward*by, c.reverse*by)

type
  # Store alignment ends by value so the hot coverage loop does not allocate
  # one heap object per read pushed into the priority queue.
  covEnd = object
    stop: pos_t
    reverse: bool

proc topStop(q: HeapQueue[covEnd]): pos_t {.inline.} = q[0].stop
proc topReverse(q: HeapQueue[covEnd]): bool {.inline.} = q[0].reverse
proc empty(q: HeapQueue[covEnd]): bool {.inline.} = len(q) == 0

proc `<`(a, b: covEnd): bool = a.stop < b.stop

proc `$`(c: coverage_t): string =
  "c=" & $(c.forward + c.reverse) & "(" & $c.forward & "+/" & $c.reverse & "-)"

type
  coverage_interval_t = genomic_interval_t[tuple[l1: coverage_t, l2: string]] # l2 is the target interval or the chromosome

proc getReadEnd(start, stop, reflen: pos_t, extend: int): pos_t =
  if extend == 0:
    return stop
  elif start < stop:
    if start + extend > reflen:
      return reflen
    else:
      return start + extend
  else:
    # start > stop
    if stop - extend < 0:
      return 0
    else:
      return stop - extend

proc coverage_iter(bam: Bam, opts: input_option_t, target: target_t): iterator(): coverage_interval_t =
  result = iterator(): coverage_interval_t {.closure.} =
    var
      next_alignment                                     = alignment_stream(bam, opts, target)
      next_change             : pos_t                    = 0
      aln                     : genomic_interval_t[bool] = next_alignment()
      more_alignments         : bool                     = not finished(next_alignment)
      target_idx              : target_index_t

    for reference in bam.hdr.targets():
      let
        reflen: pos_t = pos_t(reference.length)
        refname = reference.name
        refid: chrom_t = reference.tid

      dbEcho("new reference start:", refname, ",", $reflen, "bp")

      var
        last_pos : pos_t = 0
        coverage_ends = initHeapQueue[covEnd]()
        cov           = newCov()
        more_alignments_for_ref = more_alignments and aln.chrom == refid

      while true:
        dbEcho("Aln:", if more_alignments: $aln else: "no more", "in", refname)

        # calculate the position of the next coverage change
        if more_alignments_for_ref:
          if coverage_ends.empty():
            next_change = aln.start
          else:
            next_change = min(aln.start, coverage_ends.topStop())
        else:
          if coverage_ends.empty():
            next_change = reflen
          else:
            next_change = coverage_ends.topStop()

        if debugCov: stderr.writeLine("Last pos: " & $last_pos & ", next pos: " & $next_change)
        covAssert(last_pos < next_change or (last_pos == 0 and next_change == 0),
          "coverage went backwards from " & $last_pos & " to " & $next_change & ", at " & refname & ":" & $aln.start)
        covAssert(coverage_ends.len() == cov.tot(), "coverage not equal to queue size")
        if len(target) == 0:
          yield (refid, last_pos, next_change, (cov, refname))
        else:
          dbEcho("pre-intersection coverage:", $(refname, last_pos, next_change, cov))
          for i in intersections((refid, last_pos, next_change, cov), target, target_idx):
            yield i

        if next_change == reflen:
          break

        if debugCov:
          stderr.writeLine("-+-  next=", next_change, "\tMoreAln=", more_alignments, "|", more_alignments_for_ref, ";Cov=", cov.tot(), ";Size=", len(coverage_ends))
          if more_alignments:
            stderr.writeLine(" +-> more aln @ chr=", aln.chrom, ",pos=", aln.start)

        while more_alignments_for_ref and (next_change == aln.start):
          let readEnd = getReadEnd(aln.start, aln.stop, reflen, opts.extendFrag)

          coverage_ends.push(covEnd(stop: readEnd, reverse: aln.label))
          cov.inc(aln.label)
          if debugCov: stderr.writeLine("Added aln: " & $aln)

          aln = next_alignment()
          more_alignments = not finished(next_alignment)
          more_alignments_for_ref = more_alignments and aln.chrom == refid

        while not coverage_ends.empty() and next_change == coverage_ends.topStop():
          cov.dec(coverage_ends.topReverse())
          discard coverage_ends.pop()

        if last_pos == reflen:
          covAssert(cov.tot() == 0, "coverage not null at the end of chromosome " & refname & ": cov.tot=" & $cov.tot() & " = For:" & $cov.forward & "+Rev:" & $cov.reverse)
          covAssert(coverage_ends.len() == 0, "coverage queue not null at the end of chromosome " & refname & ": " & $coverage_ends.len())
          break

        last_pos = next_change

      if coverage_ends.len() != 0:
        stderr.writeLine("ERROR: Coverage not zero when expected. Try samtools fixmate.")
        quit(1)
    covAssert(not more_alignments, "Is the BAM sorted?")

type
  output_format_t = enum # supported output formats
    of_bed,              # BedGraph-like interval coverage format
    of_wig_fixstep,      # fixed step WIG format
  span_func_t = enum sf_max, sf_min, sf_mean # how to summarize coverage in a WIG span
  output_option_t = tuple[
    strand: bool,         # output strand-specific coverage and stats
    # coverage specific options
    no_coverage: bool,    # do not output coverage, only stats
    quantization: string, #
    output_format: output_format_t, # output format
    output_path: string,  # output file, or "-" for stdout
    span_length: pos_t,   # span for wig output format
    span_func: span_func_t,
    # stats specific options
    low_cov: int64        # report length of regions under low_cov in stats
  ]

  output_t = object
    queued: genomic_interval_t[coverage_t]
    # data fields needed for spanned output (wig)
    current_span: genomic_interval_t[coverage_t]

    opts: output_option_t
    quantization_index2label: seq[string]
    quantization_coverage2index: seq[int]
    chrom2str: TableRef[chrom_t, string]
    chrom2len: TableRef[chrom_t, pos_t]
    outFile: File

proc newSpanValue(span_func: span_func_t): coverage_t =
  case span_func:
    of sf_max, sf_mean: newCov(0, 0)
    of sf_min: newCov(high(int), high(int))

proc wigSpanWidth(o: output_t, span: genomic_interval_t[coverage_t]): pos_t =
  # The last fixedStep bin can extend past the contig end; clip it before
  # computing mean coverage so the partial bin is normalized correctly.
  min(span.stop, o.chrom2len[span.chrom]) - span.start

proc output_wig_span(o: output_t, span: genomic_interval_t[coverage_t]) =
  let span_length = float(o.wigSpanWidth(span))
  let value_str =
    if o.opts.strand:
      case o.opts.span_func:
        of sf_max, sf_min: $span.label.forward & "\t" & $span.label.reverse
        of sf_mean:
          let mean = span.label/float(span_length)
          $mean.forward & "\t" & $mean.reverse
    else:
      let tot = span.label.forward + span.label.reverse
      case o.opts.span_func:
        of sf_max, sf_min: $tot
        of sf_mean: $(float(tot)/float(span_length))
  o.outFile.writeLine(value_str)

proc write_output(o: var output_t, i: genomic_interval_t[coverage_t]) =
  if o.current_span.chrom != i.chrom or i.start < i.stop: # skip empty intervals
    case o.opts.output_format:
      of of_bed: # bed format output
        let interval_str = o.chrom2str[i.chrom] & "\t" & $i.start & "\t" & $i.stop & "\t"
        let coverage_str =
          if o.opts.strand:
            if len(o.quantization_index2label) > 0:
              o.quantization_index2label[i.label.forward] & "\t" & o.quantization_index2label[i.label.reverse]
            else:
              $i.label.forward & "\t" & $i.label.reverse
          else:
            if len(o.quantization_index2label) > 0:
              o.quantization_index2label[i.label.forward]
            else:
              $i.label.forward
        o.outFile.writeLine(interval_str & coverage_str)
      of of_wig_fixstep:
        if len(o.quantization_index2label) > 0:
          stderr.writeLine("ERROR: wig output does not support quantized coverage")
          quit(1)
        let span_length = o.opts.span_length
        if o.current_span.chrom != i.chrom: # start new contig
          if o.current_span.chrom != -1 and o.current_span.start < o.chrom2len[o.current_span.chrom]: # output last possibly incomplete span from previous chrom
            o.output_wig_span(o.current_span)
          if i.chrom == -1:
            return

          o.current_span.chrom = i.chrom
          o.current_span.start = 0
          o.current_span.stop = o.current_span.start + span_length
          # Reset the full accumulated span state on every contig switch so
          # stranded WIG output cannot inherit reverse coverage from the
          # previous chromosome.
          o.current_span.label = newSpanValue(o.opts.span_func)
          o.outFile.writeLine("fixedStep chrom=" & o.chrom2str[o.current_span.chrom] & " start=1 step=" & $span_length & " span=" & $span_length)

        while o.current_span.start <= i.stop:
          let inter = intersection_first(o.current_span, i)
          if not is_empty(inter): # update the current span value
            o.current_span.label = case o.opts.span_func:
              of sf_max: max(o.current_span.label, i.label)
              of sf_min: min(o.current_span.label, i.label)
              of sf_mean: o.current_span.label + i.label*int(len(inter))
          if inter.stop == o.current_span.stop: # span is concluded
            o.output_wig_span(o.current_span)
            # next span
            o.current_span.start += span_length
            o.current_span.stop = o.current_span.start + span_length
            o.current_span.label = newSpanValue(o.opts.span_func)
          else: # span extends beyond the interval, we are done
            break

proc push_interval(o: var output_t, i: coverage_interval_t) =
  let q = o.queued
  var c: coverage_t = i.label.l1

  dbEcho("push_interval: ", $i)
  # handle stranded output
  if not o.opts.strand:
    c.forward = c.forward + c.reverse
    c.reverse = 0
  # handle quantized output
  let qmax_cov = len(o.quantization_coverage2index)
  if qmax_cov > 0:
    c.forward = o.quantization_coverage2index[min(c.forward, qmax_cov - 1)]
    if o.opts.strand:
      c.reverse = o.quantization_coverage2index[min(c.reverse, qmax_cov - 1)]

  if q.label == c and q.stop == i.start and q.chrom == i.chrom: # extend previous interval
    if debugCov: stderr.writeLine("push_inteval: extend " & $q)
    o.queued.stop = i.stop
  else: # output previous interval and queue the new one
    o.write_output(q)
    o.queued = (i.chrom, i.start, i.stop, c)

proc flush_output(o: var output_t) =
  # Flush explicitly at the end of bam2stats instead of relying on a
  # destructor with side effects.
  case o.opts.output_format:
    of of_bed: o.write_output(o.queued)
    of of_wig_fixstep:
      o.write_output(o.queued)
      o.write_output((chrom_t(-1), pos_t(0), pos_t(0), newCov()))
  if o.opts.output_path != "-":
    o.outFile.close()

# output quantization
proc parse_quantization(o: var output_t, breaks: string) =
  var bb: seq[int] = @[0]
  for b in split(breaks, ','):
    bb.add(parse_int(b))
  bb = deduplicate(sorted(bb), isSorted = true)
  let max_right = max(bb)
  for i in 0..(len(bb) - 2):
    let
      left = bb[i]
      right = bb[i + 1] - 1
    o.quantization_index2label.add($left & "-" & $right)
    for c in left..right:
      o.quantization_coverage2index.add(i)
  o.quantization_index2label.add($max_right & "-")
  o.quantization_coverage2index.add(len(o.quantization_index2label) - 1)

proc newOutput(opts: output_option_t, bam: Bam): output_t =
  var o = output_t(
    opts: opts,
    queued: (chrom_t(-1), pos_t(0), pos_t(0), newCov()),
    quantization_index2label: @[],
    quantization_coverage2index: @[],
    chrom2str: newTable[chrom_t, string](),
    chrom2len: newTable[chrom_t, pos_t](),
    current_span: (chrom_t(-1), pos_t(0), pos_t(0), newCov()),
    outFile: if opts.no_coverage or opts.output_path == "-": stdout else: open(opts.output_path, fmWrite)
  )
  for t in bam.hdr.targets:
    o.chrom2str[t.tid] = t.name
    o.chrom2len[t.tid] = pos_t(t.length)
  if opts.quantization != "nil":
    o.parse_quantization(opts.quantization)
  o

type
  coverage_stats_t[T] = tuple[total, forward, reverse: T]
  interval_stats_t = tuple[bases, min_cov, max_cov: coverage_stats_t[int64], low_length: coverage_stats_t[pos_t], length: pos_t]
  target_stat_t = ref object
    opts: output_option_t
    stats: TableRef[string, interval_stats_t]

proc `max`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (max(s1.total, s2.total), max(s1.forward, s2.forward), max(s1.reverse, s2.reverse))
proc `min`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (min(s1.total, s2.total), min(s1.forward, s2.forward), min(s1.reverse, s2.reverse))
proc `+`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (s1.total + s2.total, s1.forward + s2.forward, s1.reverse + s2.reverse)
proc `+`[T](s: coverage_stats_t[T], x: T): coverage_stats_t[T] =
  (s.total + x, s.forward + x, s.reverse + x)
proc `*`[T](s: coverage_stats_t[T], x: T): coverage_stats_t[T] =
  (s.total*x, s.forward*x, s.reverse*x)

proc new_stats(opts: output_option_t): target_stat_t =
  target_stat_t(opts: opts, stats: newTable[string, interval_stats_t]())

proc push_interval(self: var target_stat_t, i: coverage_interval_t) =
  # update coverage statistics
  let
    l = len(i)
    name: string = i.label.l2 # target interval name
    cov: coverage_stats_t[int64] = (int64(i.label.l1.forward + i.label.l1.reverse), int64(i.label.l1.forward), int64(i.label.l1.reverse))
    low_cov = self.opts.low_cov
    low_length: coverage_stats_t[int64] = (
      (if cov.total   < low_cov: int64(l) else: 0),
      (if cov.forward < low_cov: int64(l) else: 0),
      (if cov.reverse < low_cov: int64(l) else: 0)
    )
  self.stats[name] =
    if name in self.stats:
      let o = self.stats[name]
      (
        bases: o.bases + (cov*l),
        min_cov: min(o.min_cov, cov),
        max_cov: max(o.max_cov, cov),
        low_length: o.low_length + low_length,
        length: o.length + l
      )
    else:
      (
        bases: cov*l,
        min_cov: cov,
        max_cov: cov,
        low_length: low_length,
        length: l
      )

proc mean(s: interval_stats_t): coverage_stats_t[float] =
  let l = float(s.length)
  (float(s.bases.total)/l, float(s.bases.forward)/l, float(s.bases.reverse)/l)

proc to_string[T](s: coverage_stats_t[T], strand: bool = true, sep: string = " "): string =
  if strand:
    $s.total & sep & $s.forward & sep & $s.reverse
  else:
    $s.total

proc stat_columns(self: target_stat_t, sep: string = " ", prefix: string = ""): string =
  let
    cov_cols = @["bases", "mean", "min", "max"] & (if self.opts.low_cov > 0: @["low" & $self.opts.low_cov] else: @[])
    strand_suffixes = if self.opts.strand: @["", "_forward", "_reverse"] else: @[""]
  var r: string = ""
  for mid in cov_cols:
    for suf in strand_suffixes:
      r = r & sep & prefix & mid & suf
  r & sep & prefix & "length"

proc to_string(self: target_stat_t, name: string, sep: string = " "): string =
  let strand = self.opts.strand
  if name in self.stats:
    let s = self.stats[name]
    var r =
      to_string(s.bases, strand, sep) & sep &
      to_string(s.mean, strand, sep) & sep &
      to_string(s.min_cov, strand, sep) & sep &
      to_string(s.max_cov, strand, sep) & sep
    if self.opts.low_cov > 0:
      r = r & to_string(s.low_length, strand, sep) & sep
    r & $s.length
  else:
    # This is reached when a target interval received no alignments.
    let nullstring = "0\t"
    var r = nullstring.repeat((1 + (if self.opts.low_cov > 0: 5 else: 4)*(if strand: 3 else: 1)))
    # remove last tab from string
    r[0 .. ^2]

# process coverage from a single file:
# open bam, compute coverage, print coverage output (based on outopts) and return coverage stats
proc bam2stats(bam_path: string, inopts: input_option_t, outopts: output_option_t, bam_threads: int = 0): tuple[stats: target_stat_t, target: target_t] =
  var
    bam: Bam
    target_stats: target_stat_t = new_stats(outopts)

  if bam_path == "-":
    stderr.writeLine("Reading from STDIN [press Ctrl-C to quit]")

  if not open(bam, bam_path, threads = bam_threads):
    stderr.writeLine("ERROR: Failed to open BAM/CRAM file: ", bam_path)
    quit(1)

  if bam.hdr.isNil:
    stderr.writeLine("ERROR: Invalid or empty BAM/CRAM file")
    quit(1)

  var target: target_t
  try:
    target = cookTarget(inopts.target, bam)
  except ValueError as e:
    stderr.writeLine("ERROR: ", e.msg)
    quit(1)
  var output: output_t = newOutput(outopts, bam)

  var cov_iter = coverage_iter(bam, inopts, target)
  for cov_inter in cov_iter():
    dbEcho("coverage:", $cov_inter)
    target_stats.push_interval(cov_inter)
    if not outopts.no_coverage:
      output.push_interval(cov_inter)
  if not outopts.no_coverage:
    output.flush_output()
  (target_stats, target)

proc cov*(argv: var seq[string]): int =
  let doc = """
BamToCov $version - per-base / per-target coverage

Usage:
  bamto cov [options] [<bam>...]

Input:
  <bam>...                    BAM/CRAM/SAM input files. Use - or omit for STDIN.
  -T, --threads <N>           BAM decompression threads [default: 0]
  -Q, --min-mapq <N>          Minimum mapping quality [default: 0]
  -F, --exclude-flag <FLAG>   Exclude reads with any SAM flag bit set [default: 1796]

Coverage:
  -s, --stranded              Report forward and reverse coverage separately
  -p, --physical              Calculate fragment/physical coverage
  --extend-reads <BP>         Extend reads to BP bases [default: 0]

Targets:
  -r, --regions <FILE>        BED/GFF/GTF target intervals
  --target-format <FMT>       auto|bed|gff|gtf [default: auto]
  --gff-feature <TYPE>        GFF/GTF feature type [default: CDS]
  --gff-attribute <KEY>       GFF/GTF attribute used as interval name [default: ID]
  --gff-separator <SEP>       GFF/GTF attribute separator [default: ;]

Output:
  -f, --format <FMT>          bedgraph|wig [default: bedgraph]
  -o, --output <FILE>         Coverage output file [default: -]
  --summary <FILE>            Write per-contig/target summary TSV
  --summary-only              Write summary without coverage output
  --low-cov <N>               Report bases with coverage below N [default: 0]
  --quantize <BREAKS>         Quantize coverage using comma-separated breaks

WIG options:
  --span <BP>                 WIG fixed-step span [default: 1]
  --span-op <FUNC>            mean|min|max [default: max]

Other:
  --debug                     Enable diagnostics
  -h, --help                  Show this help
""".replace("$version", version())

  let removedOptions = {
    "--report": "--summary",
    "--skip-output": "--summary-only",
    "--wig": "--format wig --span",
    "-w": "--format wig --span",
    "--op": "--span-op",
    "--mapq": "--min-mapq",
    "--flag": "--exclude-flag",
    "--gff-type": "--gff-feature",
    "-t": "--gff-feature",
    "--gff-id": "--gff-attribute",
    "-i": "--gff-attribute",
    "--gff": "--target-format gff",
    "--gtf": "--target-format gtf",
    "--report-low": "--low-cov",
    "-q": "--quantize",
  }.toTable

  for arg in argv:
    let opt = if "=" in arg: arg.split("=", 1)[0] else: arg
    if opt in removedOptions:
      stderr.writeLine("ERROR: ", opt, " was removed from 'bamto cov'; use ", removedOptions[opt], ".")
      return 1

  var normalizedArgv = argv
  for arg in normalizedArgv.mitems:
    if arg == "--extendReads":
      arg = "--extend-reads"
    elif arg.startsWith("--extendReads="):
      arg = "--extend-reads=" & arg.substr("--extendReads=".len)

  let args = docopt(doc, version = version(), argv = @["cov"] & normalizedArgv)

  debugCov = bool(args["--debug"])
  if debugCov:
    dbEcho("args:", $args)

  let
    threads = parse_int($args["--threads"])
    target_file = $args["--regions"]
    target_format = $args["--target-format"]
    coverage_format = $args["--format"]
    span_op = $args["--span-op"]
    span_length = parse_int($args["--span"])

  var
    format_gff = false
    format_gtf = false

  if target_format notin @["auto", "bed", "gff", "gtf"]:
    stderr.writeLine("ERROR: --target-format must be one of auto, bed, gff, gtf; got: ", target_format)
    quit(1)

  if coverage_format notin @["bedgraph", "wig"]:
    stderr.writeLine("ERROR: --format must be one of bedgraph, wig; got: ", coverage_format)
    quit(1)

  if span_op notin @["mean", "min", "max"]:
    stderr.writeLine("ERROR: --span-op must be one of mean, min, max; got: ", span_op)
    quit(1)

  if span_length <= 0:
    stderr.writeLine("ERROR: --span must be greater than 0; got: ", span_length)
    quit(1)

  if target_file != "nil":
    if not fileExists(target_file):
      stderr.writeLine("ERROR: Target file not found:", target_file)
      quit(1)

    case target_format:
      of "gff":
        format_gff = true
      of "gtf":
        format_gtf = true
      of "auto":
        if target_file.toLower().contains(".gff"):
          dbEcho("Parsing target as GFF")
          format_gff = true
        elif target_file.toLower().contains(".gtf"):
          dbEcho("Parsing target as GTF")
          format_gtf = true
        else:
          dbEcho("Parsing target as BED")
      else:
        dbEcho("Parsing target as BED")

  let
    gffField = $args["--gff-feature"]
    gffSeparator = $args["--gff-separator"]
    gffIdentifier = $args["--gff-attribute"]
    target_order =
      if target_file != "nil":
        target_names_in_order(target_file, format_gff, format_gtf, gffField, gffSeparator, gffIdentifier)
      else:
        @[]

  let
    input_paths =
      if len(@(args["<bam>"])) > 0:
        @(args["<bam>"])
      else:
        @["-"]
    input_opts: input_option_t = (
      min_mapping_quality: uint8(parse_int($args["--min-mapq"])),
      eflag: uint16(parse_int($args["--exclude-flag"])),
      physical: bool(args["--physical"]),
      extendFrag: parse_int($args["--extend-reads"]),
      target: if format_gff: gff_to_table(target_file, gffField, gffSeparator, gffIdentifier)
              elif format_gtf: gtf_to_table(target_file, gffField, gffSeparator, gffIdentifier)
              else: bed_to_table(target_file)
    )
    output_opts: output_option_t = (
      strand: bool(args["--stranded"]),
      no_coverage: bool(args["--summary-only"]) or len(input_paths) > 1,
      quantization: $args["--quantize"],
      output_format: if coverage_format == "wig": of_wig_fixstep else: of_bed,
      output_path: $args["--output"],
      span_length: pos_t(span_length),
      span_func: if span_op == "max": sf_max
                 elif span_op == "min": sf_min
                 else: sf_mean,
      low_cov: int64(parse_int($args["--low-cov"]))
    )

  # Preflight check input files
  var missing_files = 0
  for inputBam in input_paths:
    if inputBam == "-":
      dbEcho("Will read STDIN")
    elif not fileExists(inputBam):
      missing_files += 1
      stderr.writeLine("ERROR: Input BAM file <", inputBam, "> not found.")
  if missing_files > 0:
    quit(1)

  # Multiple BAMs
  if len(input_paths) > 1:
    if not bool(args["--summary-only"]):
      stderr.writeLine("WARNING: coverage output for multiple input files is not implemented, so it will not be produced; use --summary-only to suppress this warning")
    if target_file == "nil":
      stderr.writeLine("ERROR: Multiple BAMs are handled via target file (--regions). Supply a target.")
      quit(1)

  # All the number crunching lives in bam2stats, ready to be threaded.
  var bam_stats: seq[tuple[stats: target_stat_t, target: target_t]]
  for p in input_paths:
    dbEcho("running", p)
    bam_stats.add(bam2stats(p, input_opts, output_opts, bam_threads = threads))

  if args["--summary"]: # print summary table
    dbEcho("stats reporting")
    # assemble table index
    var index = initOrderedSet[string]()
    let
      sample_names = input_paths
      target = bam_stats[0].target # get cooked target from the first bam

    if len(input_opts.target) > 0:
      # get interval names from target in the order they appear
      for name in target_order:
        index.incl(name)
    else:
      # if there is no target, use chromosomes
      for s in bam_stats:
        for t in s.stats.stats.keys():
          index.incl(t)

    # check that each interval in stats has been put in index
    for s in bam_stats:
      for t in s.stats.stats.keys():
        if not (t in index):
          covAssert(t in index, "t not in index: " & t)
    dbEcho("target:", $target)
    dbEcho("index:", $index)

    # print header
    dbEcho("report: header")
    let
      report = open($args["--summary"], fmWrite)
      sep = "\t"
    report.write("interval")
    for x in zip(sample_names, bam_stats):
      report.write(stat_columns(x[1].stats, sep = sep, prefix = x[0] & "_"))
    report.write("\n")

    # print body
    dbEcho("report: body")
    for t in index:
      report.write(t)
      for s in bam_stats:
        report.write(sep & to_string(s.stats, t, sep = sep))
      report.write("\n")

  dbEcho("exiting successfully!")
  return 0
