#!/usr/bin/env python3
"""
Generate a synthetic MinKNOW-style `bam_pass/` folder for testing Shiny-Sturgeon
without patient data or a sequencing run.

The output looks like what MinKNOW writes when it basecalls, calls
methylation and aligns in real time:

    <out>/barcode01/synthetic_0.bam, synthetic_1.bam, ...
    <out>/unclassified/synthetic_unclassified_0.bam, ...

Every read is aligned to a T2T-CHM13v2.0 contig (header only, no FASTA needed)
and carries 5mC/5hmC calls in MM/ML tags at CpG sites, so samtools and
`modkit adjust-mods` / `modkit extract full` can process the files.

The data is RANDOM. It exercises the pipeline (file discovery, merging,
modkit, sturgeon inputtobed, CNV binning) but any classification made from it
is meaningless. Never use it to judge classifier performance.

Examples
--------
    # 5 BAM files of 400 reads each, with CpG calls at real Sturgeon probe
    # positions and a fake chr7 gain / chr10 loss for the CNV plot
    python tools/make_synthetic_bams.py --out /tmp/demo/bam_pass \\
        --probes /path/to/probelocs_chm13.bed --gain 7 --loss 10

    # Simulate a live run: write one file every 30 seconds while the app
    # is monitoring the folder
    python tools/make_synthetic_bams.py --out /tmp/demo/bam_pass --live-interval 30

Needs only Python 3 and samtools (both in the `sturgeon` conda environment
that install.sh creates): reads are written as SAM text and converted to BAM
by samtools.
"""

import argparse
import bisect
import os
import random
import shutil
import subprocess
import sys
import time

# T2T-CHM13v2.0 contig lengths (chm13v2.0.fa)
CHM13_LENGTHS = {
    "1": 248387328, "2": 242696752, "3": 201105948, "4": 193574945,
    "5": 182045439, "6": 172126628, "7": 160567428, "8": 146259331,
    "9": 150617247, "10": 134758134, "11": 135127769, "12": 133324548,
    "13": 113566686, "14": 101161492, "15": 99753195, "16": 96330374,
    "17": 84276897, "18": 80542538, "19": 61707364, "20": 66210255,
    "21": 45090682, "22": 51324926, "X": 154259566, "Y": 62460029,
    "M": 16569,
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate synthetic MinKNOW-style BAM files with methylation tags.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--out", required=True, help="bam_pass folder to create or add to")
    p.add_argument("--barcode", default="01", help="barcode number (zero-padded to two digits)")
    p.add_argument("--n-files", type=int, default=5, help="BAM files to write for the barcode")
    p.add_argument("--reads-per-file", type=int, default=400, help="reads per BAM file")
    p.add_argument("--unclassified-files", type=int, default=1,
                   help="BAM files to write in unclassified/")
    p.add_argument("--read-length", type=int, default=5000, help="mean read length (bp)")
    p.add_argument("--probes", help="Sturgeon probes file (e.g. probelocs_chm13.bed). "
                                    "If given, CpG calls are placed at probe positions; "
                                    "otherwise at random positions.")
    p.add_argument("--chrom-style", choices=["chr", "plain"], default="chr",
                   help="contig names as 'chr1' (chm13v2.0.fa) or '1'")
    p.add_argument("--gain", action="append", default=[], metavar="CHROM",
                   help="chromosome with a fake copy-number gain (1.5x reads); repeatable")
    p.add_argument("--loss", action="append", default=[], metavar="CHROM",
                   help="chromosome with a fake copy-number loss (0.5x reads); repeatable")
    p.add_argument("--no-sequence", action="store_true",
                   help="omit read sequences and methylation tags (much smaller files; "
                        "enough for CNV plots, not for modkit)")
    p.add_argument("--seed", type=int, default=1, help="random seed (same seed, same files)")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--live-interval", type=float, default=0,
                   help="seconds to wait between files, to simulate a live run")
    return p.parse_args(argv)


def load_probes(path, wanted_chroms):
    """Return {plain chrom: sorted list of 0-based CpG positions} from a Sturgeon probes file."""
    probes = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 2 or not fields[1].isdigit():
                continue  # header or malformed line
            chrom = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            if chrom in wanted_chroms:
                probes.setdefault(chrom, []).append(int(fields[1]))
    for positions in probes.values():
        positions.sort()
    if not probes:
        sys.exit(f"No usable probe positions found in {path}")
    return probes


class ReadFactory:
    def __init__(self, args):
        self.args = args
        self.rng = random.Random(args.seed)
        self.chroms = [c for c in CHM13_LENGTHS if c != "M"]
        weights = []
        for c in self.chroms:
            w = CHM13_LENGTHS[c]
            if c in args.gain:
                w *= 1.5
            if c in args.loss:
                w *= 0.5
            weights.append(w)
        self.weights = weights
        self.probes = load_probes(args.probes, set(self.chroms)) if args.probes else None
        # Each CpG site gets a fixed methylation level (bimodal, like real
        # tissue) so repeated reads over the same site roughly agree.
        self.site_seed = args.seed * 1_000_003

    def contig_name(self, chrom):
        return f"chr{chrom}" if self.args.chrom_style == "chr" else chrom

    def header(self):
        lines = ["@HD\tVN:1.6\tSO:coordinate"]
        lines += [f"@SQ\tSN:{self.contig_name(c)}\tLN:{CHM13_LENGTHS[c]}" for c in CHM13_LENGTHS]
        lines.append("@PG\tID:make_synthetic_bams\tPN:make_synthetic_bams.py\tCL:" + " ".join(sys.argv))
        return "\n".join(lines) + "\n"

    def site_beta(self, chrom, pos):
        r = random.Random(f"{self.site_seed}:{chrom}:{pos}").random()
        return 0.1 if r < 0.5 else 0.85

    def cpg_offsets(self, chrom, start, length):
        """0-based offsets within the read where a CpG 'C' sits."""
        if self.probes is not None:
            positions = self.probes.get(chrom, [])
            lo = bisect.bisect_left(positions, start)
            hi = bisect.bisect_left(positions, start + length - 1)
            offsets = [p - start for p in positions[lo:hi]]
        else:
            offsets = []
            pos = self.rng.randint(20, 300)
            while pos < length - 1:
                offsets.append(pos)
                pos += self.rng.randint(50, 400)
        # keep sites at least 2 bp apart so the forced "CG" dinucleotides don't overlap
        kept, last = [], -10
        for o in offsets:
            if o - last >= 2 and o < length - 1:
                kept.append(o)
                last = o
        return kept

    def make_read(self, name):
        """Return (sort key, SAM line) for one forward-strand, fully aligned read."""
        chrom = self.rng.choices(self.chroms, weights=self.weights)[0]
        length = max(200, int(self.rng.gauss(self.args.read_length, self.args.read_length * 0.3)))
        start = self.rng.randint(0, CHM13_LENGTHS[chrom] - length - 1)
        key = (list(CHM13_LENGTHS).index(chrom), start)
        # flag 0 = forward strand, primary, mapped; SAM positions are 1-based
        fields = [name, "0", self.contig_name(chrom), str(start + 1), "60", f"{length}M", "*", "0", "0"]

        if self.args.no_sequence:
            return key, "\t".join(fields + ["*", "*"])

        offsets = self.cpg_offsets(chrom, start, length)
        # random sequence without any C, then a "CG" at every called site, so
        # the only Cs in the read are the called ones and the MM skip counts are 0
        seq = [self.rng.choice("AGT") for _ in range(length)]
        for o in offsets:
            seq[o], seq[o + 1] = "C", "G"
        fields += ["".join(seq), "5" * length]  # quality "5" = Q20

        if offsets:
            h_probs, m_probs = [], []
            # h and m are probabilities of the same C (scaled to 0-255), so
            # h + m must stay <= 255 or modkit rejects the read
            for o in offsets:
                methylated = self.rng.random() < self.site_beta(chrom, start + o)
                if methylated and self.rng.random() < 0.1:  # a few 5hmC calls for adjust-mods
                    h = self.rng.randint(180, 230)
                    m = self.rng.randint(0, min(40, 255 - h))
                elif methylated:
                    h = self.rng.randint(0, 20)
                    m = self.rng.randint(200, 255 - h)
                else:
                    h, m = self.rng.randint(0, 20), self.rng.randint(0, 30)
                h_probs.append(h)
                m_probs.append(m)
            skips = ",".join(["0"] * len(offsets))
            fields.append(f"MM:Z:C+h?,{skips};C+m?,{skips};")
            fields.append("ML:B:C," + ",".join(map(str, h_probs + m_probs)))  # uint8 probabilities
        return key, "\t".join(fields)


def write_bam(factory, path, n_reads, read_prefix):
    reads = [factory.make_read(f"{read_prefix}_{i:06d}") for i in range(n_reads)]
    reads.sort(key=lambda r: r[0])
    sam = factory.header() + "".join(line + "\n" for _, line in reads)
    # write under a hidden name first (R's list.files() skips dotfiles), so
    # the app never picks up a half-written file
    tmp = os.path.join(os.path.dirname(path), "." + os.path.basename(path) + ".tmp")
    subprocess.run([factory.args.samtools, "view", "-b", "-o", tmp, "-"], input=sam.encode(), check=True)
    os.replace(tmp, path)  # appear in the folder only once complete, like a finished MinKNOW file
    subprocess.run([factory.args.samtools, "index", path], check=True)


def main(argv=None):
    args = parse_args(argv)
    if shutil.which(args.samtools) is None:
        sys.exit(f"samtools not found ({args.samtools}); pass its full path with --samtools")
    args.barcode = args.barcode.zfill(2)
    args.gain = [c.replace("chr", "") for c in args.gain]
    args.loss = [c.replace("chr", "") for c in args.loss]
    factory = ReadFactory(args)

    barcode_dir = os.path.join(args.out, f"barcode{args.barcode}")
    unclassified_dir = os.path.join(args.out, "unclassified")
    os.makedirs(barcode_dir, exist_ok=True)
    os.makedirs(unclassified_dir, exist_ok=True)

    jobs = [(barcode_dir, f"synthetic_{i}") for i in range(args.n_files)]
    jobs += [(unclassified_dir, f"synthetic_unclassified_{i}") for i in range(args.unclassified_files)]

    for n, (folder, stem) in enumerate(jobs):
        if n > 0 and args.live_interval > 0:
            time.sleep(args.live_interval)
        path = os.path.join(folder, stem + ".bam")
        write_bam(factory, path, args.reads_per_file, f"{os.path.basename(folder)}_{stem}")
        print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
