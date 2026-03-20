import std/[strutils, tables, math, algorithm]
import lapper
import hts/files
import cligen


type
  ExonIv* = object
    startPos*, stopPos*: int
    tx_id*: string

  TxStats = object
    tin, minCov, maxCov, meanCov, medianCov: float
    length: int


proc start*(iv: ExonIv): int {.inline.} = iv.startPos
proc stop*(iv: ExonIv): int {.inline.} = iv.stopPos


proc getTranscriptId(attrStr: string): string =
  for item in attrStr.split(';'):
    let pair = item.strip().split(' ')
    if pair.len == 2 and pair[0] == "transcript_id":
      return pair[1].replace("\"", "")
  return ""


proc buildGtfIndex(gtfFile: string): tuple[trees: Table[string, Lapper[ExonIv]], txLengths: Table[string, int]] =
  var 
    exonDict = initTable[string, seq[ExonIv]]()
    txLengths = initTable[string, int]()

  for line in lines(gtfFile):
    if line.startsWith("#"): continue
    let cols = line.split('\t')
    if cols.len < 9 or cols[2] != "exon": continue

    let chrom = cols[0]
    let start = parseInt(cols[3]) - 1
    let stop = parseInt(cols[4])
    let txId = getTranscriptId(cols[8])
    
    if txId.len > 0:
      let length = stop - start
      txLengths[txId] = txLengths.getOrDefault(txId, 0) + length
      
      if not exonDict.hasKey(chrom):
        exonDict[chrom] = @[]
      exonDict[chrom].add(ExonIv(startPos: start, stopPos: stop, tx_id: txId))

  var trees = initTable[string, Lapper[ExonIv]]()
  for chrom, exons in exonDict.mpairs:
    trees[chrom] = lapify(exons)
    
  return (trees, txLengths)


proc parseBedCoverage(bedFile: string, trees: var Table[string, Lapper[ExonIv]]): tuple[totalCov: Table[string, float], logSum: Table[string, float], covBlocks: Table[string, seq[tuple[c: float, w: int]]]] =
  var 
    txTotalCov = initTable[string, float]()
    txLogSum = initTable[string, float]()
    txCovBlocks = initTable[string, seq[tuple[c: float, w: int]]]()
    f: HTSFile

  try:
    discard f.open(bedFile, "r")
  except Exception:
    quit("Error: Could not open BED file " & bedFile)

  var line = newStringOfCap(2048)
  while f.readLine(line):
    let cols = line.split('\t')
    if cols.len < 4: continue
    
    let chrom = cols[0]
    let bStart = parseInt(cols[1])
    let bStop = parseInt(cols[2])
    let cov = parseFloat(cols[3])

    if cov == 0.0 or not trees.hasKey(chrom): 
      continue
      
    for exon in trees[chrom].find(bStart, bStop):
      let txId = exon.tx_id
      let overlapStart = max(bStart, exon.startPos)
      let overlapStop = min(bStop, exon.stopPos)
      let width = float(overlapStop - overlapStart)

      if width > 0.0:
        txTotalCov[txId] = txTotalCov.getOrDefault(txId, 0.0) + (width * cov)
        txLogSum[txId] = txLogSum.getOrDefault(txId, 0.0) + (width * cov * ln(cov))
        
        # store the RLE block 
        if not txCovBlocks.hasKey(txId): txCovBlocks[txId] = @[]
        txCovBlocks[txId].add((cov, int(width)))

  f.close()
  return (txTotalCov, txLogSum, txCovBlocks)


proc computeTinScores(txLengths: Table[string, int], txTotalCov: Table[string, float], txLogSum: Table[string, float], txCovBlocks: Table[string, seq[tuple[c: float, w: int]]], minCov: float): Table[string, TxStats] =
  var results = initTable[string, TxStats]()
  
  for txId, L in txLengths:
    let C = txTotalCov.getOrDefault(txId, 0.0)
    let logSum = txLogSum.getOrDefault(txId, 0.0)
    let blocks = txCovBlocks.getOrDefault(txId, @[])
    
    var stat = TxStats(length: L, tin: 0.0, minCov: 0.0, maxCov: 0.0, meanCov: 0.0, medianCov: 0.0)
    
    if L > 0:
      stat.meanCov = C / float(L)
      
      var fullBlocks = blocks
      var covered = 0
      for b in blocks: covered += b.w
      if covered < L:
        fullBlocks.add((0.0, L - covered))
        
      fullBlocks.sort(proc(x, y: tuple[c: float, w: int]): int = cmp(x.c, y.c))
      
      if fullBlocks.len > 0:
        stat.minCov = fullBlocks[0].c
        stat.maxCov = fullBlocks[^1].c
        
        # Find the median by walking the sorted widths
        let mid1 = (L - 1) div 2
        let mid2 = L div 2
        var current = 0
        var med1, med2 = -1.0
        
        for b in fullBlocks:
          if med1 < 0.0 and current + b.w > mid1: med1 = b.c
          if med2 < 0.0 and current + b.w > mid2: med2 = b.c
          if med1 >= 0.0 and med2 >= 0.0: break
          current += b.w
          
        stat.medianCov = (med1 + med2) / 2.0

      if C >= minCov:
        var entropy = ln(C) - (logSum / C)
        if entropy < 0.0: entropy = 0.0 
        stat.tin = 100.0 * exp(entropy) / float(L)

    results[txId] = stat
    
  return results


proc main(bed: string, gtf: string, output: string = "", minCov: float = 1.0) =
  var (trees, txLengths) = buildGtfIndex(gtf)
  let (txTotalCov, txLogSum, txCovBlocks) = parseBedCoverage(bed, trees)
  let stats = computeTinScores(txLengths, txTotalCov, txLogSum, txCovBlocks, minCov)

  var outStream: File
  if output.len > 0:
    if not open(outStream, output, fmWrite):
      quit("Error: Could not open output file " & output & " for writing.")
  else:
    outStream = stdout 

  outStream.writeLine("transcript_id\tlength\tmin_cov\tmax_cov\tmean_cov\tmedian_cov\tTIN")
  for txId, s in stats:
    outStream.writeLine(txId, "\t", 
                        s.length, "\t", 
                        formatFloat(s.minCov, ffDecimal, 2), "\t",
                        formatFloat(s.maxCov, ffDecimal, 2), "\t",
                        formatFloat(s.meanCov, ffDecimal, 2), "\t",
                        formatFloat(s.medianCov, ffDecimal, 2), "\t",
                        formatFloat(s.tin, ffDecimal, 2))

  if output.len > 0:
    outStream.close()


when isMainModule:
  dispatch main, help = {
    "bed": "Path to the mosdepth per-base bed.gz file",
    "gtf": "Path to the input GTF annotations file",
    "output": "Optional path to save the TSV output. Defaults to stdout.",
    "minCov": "Minimum total accumulated depth across the transcript to calculate TIN (default: 1.0)"
  }