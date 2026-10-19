#!/usr/bin/env python3
"""Generates the BGZF / BAM / BAI / CSI / CRAM / CRAI / tabix fixtures in this directory.

The reads, variants and features are simulated here (seeded, so the output is reproducible), but every binary
file and every expected result is produced by the reference implementation, htslib/samtools 1.24, never by
Biotech.Net:

    samtools sort / index / view / idxstats   ->  reads.bam, .bai, .csi, *.cram, .crai, expected SAM and regions
    bgzip / tabix                             ->  *.gz, .tbi, .csi and expected tabix query output

Requirements: python 3.9+, and samtools, bgzip and tabix (htslib 1.24) on PATH. On Windows the MSYS2 packages
mingw-w64-ucrt-x86_64-samtools / -htslib work. Optionally pass --hts-specs <dir> (a checkout of
https://github.com/samtools/hts-specs) to refresh the conformance files under ../hts-specs.

    python gen_hts_fixtures.py [--hts-specs path/to/hts-specs]
"""
import argparse
import gzip
import hashlib
import os
import random
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = [("chr1", 40000), ("chr2", 25000), ("chr3", 8000)]  # chr3 carries no reads (empty index entry)


# CRAM: the reference is found through the @SQ UR:ref.fa tags (relative to this directory). Passing -T would make
# htslib write the absolute local path into every UR tag; REF_PATH/REF_CACHE are pointed nowhere so htslib never
# downloads a reference from the EBI server instead.
ENV = dict(os.environ, REF_PATH=os.path.join(HERE, "no-such-dir", "%s"), REF_CACHE=os.path.join(HERE, "no-such-dir", "%s"))


def run(*args, stdout=None):
    print("+", " ".join(args), file=sys.stderr)
    subprocess.run(args, check=True, cwd=HERE, stdout=stdout, env=ENV)


def capture(*args):
    return subprocess.run(args, check=True, cwd=HERE, capture_output=True, env=ENV).stdout


# --------------------------------------------------------------------------------------------- simulation
rng = random.Random(20261018)


def random_dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_reference():
    seqs = {}
    for name, length in REFS:
        s = list(random_dna(length))
        if name == "chr1":
            s[10000:10100] = "N" * 100  # an assembly gap
        seqs[name] = "".join(s)
    with open(os.path.join(HERE, "ref.fa"), "w", newline="\n") as f:
        for name, _ in REFS:
            f.write(f">{name} simulated\n")
            s = seqs[name]
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    return seqs


QBINS = [2, 12, 23, 37]


def quals(n, binned=True):
    if binned:
        return "".join(chr(33 + rng.choice(QBINS)) for _ in range(n))
    return "".join(chr(33 + rng.randint(2, 41)) for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def simulate_alignment(ref, start, read_len, allow_indel=True, clip=True):
    """Builds (read, cigar ops, ref span) by walking the reference from `start` with random edits."""
    ops = []  # list of (op, length)
    read = []
    r = start

    def push(op, n):
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + n)
        else:
            ops.append((op, n))

    left_clip = rng.randint(3, 12) if clip and rng.random() < 0.12 else 0
    right_clip = rng.randint(3, 12) if clip and rng.random() < 0.12 else 0
    if left_clip:
        read.append(random_dna(left_clip))
        push("S", left_clip)
    body = read_len - left_clip - right_clip
    indel_at = rng.randint(15, body - 15) if allow_indel and body > 40 and rng.random() < 0.25 else -1
    q = 0
    while q < body and r < len(ref):
        if q == indel_at:
            if rng.random() < 0.5:
                n = min(rng.randint(1, 4), body - q - 1)
                read.append(random_dna(n))
                push("I", n)
                q += n
            else:
                n = rng.randint(1, 6)
                push("D", n)
                r += n
            indel_at = -1
            continue
        b = ref[r]
        x = rng.random()
        if x < 0.01:
            b = rng.choice([c for c in "ACGT" if c != b])
        elif x < 0.012:
            b = "N"
        read.append(b)
        push("M", 1)
        q += 1
        r += 1
    if right_clip:
        read.append(random_dna(right_clip))
        push("S", right_clip)
    return "".join(read), ops, r - start


def md_nm(ref, start, read, ops):
    md, nm, run_len, r, q = [], 0, 0, start, 0
    for op, n in ops:
        if op in "M=X":
            for _ in range(n):
                if read[q] == ref[r] and read[q] != "N":
                    run_len += 1
                else:
                    md.append(str(run_len))
                    md.append(ref[r])
                    run_len = 0
                    nm += 1
                q += 1
                r += 1
        elif op == "I":
            q += n
            nm += n
        elif op == "S":
            q += n
        elif op == "D":
            md.append(str(run_len))
            md.append("^" + ref[r:r + n])
            run_len = 0
            r += n
            nm += n
        elif op == "N":
            r += n
    md.append(str(run_len))
    return "".join(md), nm


def cigar_str(ops):
    return "".join(f"{n}{op}" for op, n in ops) if ops else "*"


def make_sam(seqs):
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    for name, length in REFS:
        # M5/UR given up front so htslib does not insert absolute local paths when writing CRAM
        m5 = hashlib.md5(seqs[name].upper().encode()).hexdigest()
        lines.append(f"@SQ\tSN:{name}\tLN:{length}\tM5:{m5}\tUR:ref.fa")
    lines.append("@RG\tID:grpA\tSM:sample1\tPL:ILLUMINA\tLB:lib1")
    lines.append("@RG\tID:grpB\tSM:sample1\tPL:ILLUMINA\tLB:lib2")
    lines.append("@PG\tID:sim\tPN:gen_hts_fixtures.py\tVN:1")
    lines.append("@CO\tSimulated reads for Biotech.Net BAM/CRAM tests")
    recs = []

    def rec(qname, flag, rname, pos, mapq, ops, rnext, pnext, tlen, seq, qual, tags):
        recs.append("\t".join([qname, str(flag), rname, str(pos), str(mapq), cigar_str(ops), rnext, str(pnext),
                               str(tlen), seq or "*", qual or "*"] + tags))

    weighted = ["chr1"] * 3 + ["chr2"] * 2
    for i in range(320):
        chrom = rng.choice(weighted)
        ref = seqs[chrom]
        L1 = rng.choice([100, 100, 100, 75, 150])
        L2 = rng.choice([100, 100, 75])
        frag = max(L1, L2) + int(rng.gauss(220, 40))
        start = rng.randint(0, len(ref) - frag - 10)
        rg = "grpA" if i % 3 else "grpB"
        name = f"pair{i:04d}"
        s1, o1, span1 = simulate_alignment(ref, start, L1)
        p2 = start + frag - L2
        s2, o2, span2 = simulate_alignment(ref, p2, L2)
        q1 = quals(len(s1), binned=i % 7 != 0)
        q2 = quals(len(s2), binned=i % 7 != 0)
        kind = rng.random()
        md1, nm1 = md_nm(ref, start, s1, o1)
        md2, nm2 = md_nm(ref, p2, s2, o2)
        pos1, pos2 = start + 1, p2 + 1
        tl = (p2 + span2) - start
        mq1, mq2 = rng.choice([60, 60, 60, 37, 0, 255]), rng.choice([60, 60, 23])
        extra1 = [f"AS:i:{len(s1) - 5 * nm1}", f"XS:i:{-rng.randint(0, 300)}", f"RG:Z:{rg}"]
        extra2 = [f"AS:i:{len(s2) - 5 * nm2}", f"RG:Z:{rg}"]
        if i % 11 == 0:
            extra1 += ["XA:A:" + rng.choice("ACGT"), f"XF:f:{rng.uniform(-100, 100):.4g}", "XH:H:1AE301",
                       "XZ:Z:free text with spaces", f"XB:B:c,{rng.randint(-128, 127)},-1,0,127",
                       f"XC:B:C,0,255,{rng.randint(0, 255)}", "Xs:i:0",
                       f"Xt:B:s,-32768,32767,{rng.randint(-5, 5)}", "Xu:B:S,0,65535",
                       f"Xi:B:i,-2147483648,2147483647,{rng.randint(-9, 9)}", "XI:B:I,0,4294967295",
                       "Xf:B:f,1.5,-0.25,3.14159,1e-07,123456789", "Xe:B:c",
                       f"Xl:i:{rng.choice([4294967295, 3000000000, -2147483648, 70000, -40000, -200, 300])}"]
        if kind < 0.07:  # read2 unmapped, placed at the mate's position
            rec(name, 1 + 8 + 64, chrom, pos1, mq1, o1, "=", pos1, 0, s1, q1,
                [f"NM:i:{nm1}", f"MD:Z:{md1}"] + extra1)
            rec(name, 1 + 4 + 128, chrom, pos1, 0, [], "=", pos1, 0, random_dna(L2), q2, [f"RG:Z:{rg}"])
            continue
        if kind < 0.11:  # mates on different chromosomes
            other = "chr2" if chrom == "chr1" else "chr1"
            op = rng.randint(1, REFS[0][1] - 200 if other == "chr1" else REFS[1][1] - 200)
            s2, o2, span2 = simulate_alignment(seqs[other], op - 1, L2)
            md2, nm2 = md_nm(seqs[other], op - 1, s2, o2)
            rec(name, 1 + 32 + 64, chrom, pos1, mq1, o1, other, op, 0, s1, q1, [f"NM:i:{nm1}", f"MD:Z:{md1}"] + extra1)
            rec(name, 1 + 16 + 128, other, op, mq2, o2, chrom, pos1, 0, s2, q2,
                [f"NM:i:{nm2}", f"MD:Z:{md2}"] + extra2)
            continue
        flag1 = 1 + 2 + 32 + 64
        flag2 = 1 + 2 + 16 + 128
        if rng.random() < 0.03:
            flag1 |= 0x400
            flag2 |= 0x400
        if rng.random() < 0.02:
            flag2 |= 0x200
        rec(name, flag1, chrom, pos1, mq1, o1, "=", pos2, tl, s1, q1, [f"NM:i:{nm1}", f"MD:Z:{md1}"] + extra1)
        tags2 = [f"MD:Z:{md2}", f"NM:i:{nm2}"] + extra2 if i % 5 == 0 else [f"NM:i:{nm2}", f"MD:Z:{md2}"] + extra2
        rec(name, flag2, chrom, pos2, mq2, o2, "=", pos1, -tl, s2, "*" if i % 37 == 0 else q2, tags2)
        if i % 29 == 0:  # a secondary alignment with SEQ/QUAL "*"
            sp = rng.randint(0, len(ref) - 200)
            _, so, _ = simulate_alignment(ref, sp, L1, allow_indel=False, clip=False)
            rec(name, 1 + 64 + 256, chrom, sp + 1, 0, so, "=", pos2, 0, "*", "*", ["NM:i:3", f"RG:Z:{rg}"])
        if i % 31 == 0:  # a supplementary alignment with hard clips
            sp = rng.randint(0, len(ref) - 200)
            ss, so, _ = simulate_alignment(ref, sp, 60, allow_indel=False, clip=False)
            so = [("H", L1 - 60)] + so
            mds, nms = md_nm(ref, sp, ss, so)
            rec(name, 1 + 64 + 2048, chrom, sp + 1, 17, so, "=", pos2, 0, ss, quals(len(ss)),
                [f"NM:i:{nms}", f"MD:Z:{mds}", f"SA:Z:{chrom},{pos1},+,{cigar_str(o1)},{mq1},{nm1};", f"RG:Z:{rg}"])

    # Long, spliced and otherwise unusual single-end reads.
    ref = seqs["chr1"]
    # spans ~20 kb through an N skip: lands in a high-level BAI bin and overlaps every query in between
    s_a, o_a, _ = simulate_alignment(ref, 1000, 80, allow_indel=False, clip=False)
    s_b, o_b, _ = simulate_alignment(ref, 21000, 80, allow_indel=False, clip=False)
    ops = o_a + [("N", 21000 - 1080)] + o_b
    md, nm = md_nm(ref, 1000, s_a + s_b, ops)
    rec("spliced1", 0, "chr1", 1001, 50, ops, "*", 0, 0, s_a + s_b, quals(160), [f"NM:i:{nm}", f"MD:Z:{md}", "XS:A:+"])
    for k in range(6):
        p = rng.randint(0, len(ref) - 4000)
        s, o, _ = simulate_alignment(ref, p, rng.randint(1500, 3000), clip=False)
        md, nm = md_nm(ref, p, s, o)
        rec(f"long{k}", 16 if k % 2 else 0, "chr1", p + 1, 60, o, "*", 0, 0, s, quals(len(s), binned=False),
            [f"NM:i:{nm}", f"MD:Z:{md}", "RG:Z:grpA"])
    # = / X operations, padding, a read at the very first base, lower-case bases
    s, o, _ = simulate_alignment(ref, 0, 100, allow_indel=False, clip=False)
    eqx = []
    for j, c in enumerate(s):
        eqx.append(("=" if c == ref[j] else "X", 1))
    merged = []
    for op, n in eqx:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + n)
        else:
            merged.append((op, n))
    rec("eqx_first_base", 0, "chr1", 1, 60, merged, "*", 0, 0, s, quals(100), [])
    p = 5000
    seq = ref[p:p + 40] + "AC" + "GT" + ref[p + 40:p + 80]
    rec("padded", 0, "chr1", p + 1, 60, [("M", 40), ("I", 2), ("P", 2), ("I", 2), ("M", 40)], "*", 0, 0, seq, quals(len(seq)), [])
    seq = ref[7000:7100].lower()
    rec("lowercase", 0, "chr2", 7001, 60, [("M", 100)], "*", 0, 0, seq, quals(100), [])
    # an alignment reaching the very end of chr2
    L = REFS[1][1]
    rec("end_of_chr2", 16, "chr2", L - 99, 60, [("M", 100)], "*", 0, 0, seqs["chr2"][L - 100:], quals(100), [])
    # reads over the N gap
    rec("over_gap", 0, "chr1", 9951, 3, [("M", 200)], "*", 0, 0, ref[9950:10150], quals(200), [])
    # unplaced unmapped reads (end up after all mapped reads)
    for k in range(12):
        L = rng.choice([50, 100, 151])
        tags = ["RG:Z:grpB"] if k % 2 else []
        rec(f"unplaced{k:02d}", 4, "*", 0, 0, [], "*", 0, 0, random_dna(L), quals(L) if k != 3 else "*", tags)
    rng.shuffle(recs)
    with open(os.path.join(HERE, "reads.sam"), "w", newline="\n") as f:
        f.write("\n".join(lines + recs) + "\n")


# --------------------------------------------------------------------------------------------- tabix inputs
def make_vcf(seqs):
    out = ["##fileformat=VCFv4.2"]
    for name, length in REFS:
        out.append(f"##contig=<ID={name},length={length}>")
    out += ['##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"]
    for name, length in REFS:
        pos = 0
        while True:
            pos += rng.randint(1, 40)
            if pos >= length - 60:
                break
            ref = seqs[name][pos - 1]
            k = rng.random()
            if k < 0.08:  # deletion
                n = rng.randint(2, 50)
                refa, alt = seqs[name][pos - 1:pos - 1 + n], ref
                info = f"DP={rng.randint(5, 80)}"
            elif k < 0.10:  # symbolic deletion with END
                end = min(length, pos + rng.randint(1000, 20000))
                refa, alt, info = ref, "<DEL>", f"SVTYPE=DEL;END={end}"
            else:
                refa, alt = ref, rng.choice([c for c in "ACGT" if c != ref] or ["A"])
                info = f"DP={rng.randint(5, 80)}"
            gts = "\t".join(rng.choice(["0/0", "0/1", "1/1", "./."]) for _ in range(2))
            out.append(f"{name}\t{pos}\t.\t{refa}\t{alt}\t{rng.randint(10, 99)}\tPASS\t{info}\tGT\t{gts}")
    with open(os.path.join(HERE, "variants.vcf"), "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


def make_bed():
    out = ["track name=features", "# comment line"]
    for name, length in REFS:
        rows = []
        for k in range(rng.randint(150, 300)):
            s = rng.randint(0, length - 1)
            e = min(length, s + rng.choice([1, 10, 100, 1000, 5000, 20000]))
            rows.append((s, e))
        rows.sort()
        for j, (s, e) in enumerate(rows):
            out.append(f"{name}\t{s}\t{e}\tfeat_{name}_{j}\t{rng.randint(0, 1000)}\t{rng.choice('+-')}")
    with open(os.path.join(HERE, "features.bed"), "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


def make_gff():
    out = ["##gff-version 3"]
    for name, length in REFS:
        out.append(f"##sequence-region {name} 1 {length}")
    for name, length in REFS:
        p = 1
        g = 0
        rows = []
        while p < length - 3000:
            p += rng.randint(50, 1500)
            span = rng.randint(300, 2500)
            e = min(length, p + span)
            strand = rng.choice("+-")
            gid = f"gene_{name}_{g}"
            mid = (p + e) // 2
            rows.append((p, 0, f"{name}\tsim\tgene\t{p}\t{e}\t.\t{strand}\t.\tID={gid}"))
            rows.append((p, 1, f"{name}\tsim\texon\t{p}\t{mid}\t.\t{strand}\t.\tParent={gid}"))
            rows.append((mid + 50, 2, f"{name}\tsim\texon\t{mid + 50}\t{e}\t.\t{strand}\t.\tParent={gid}"))
            g += 1
        out += [r[2] for r in sorted(rows)]
    with open(os.path.join(HERE, "genes.gff3"), "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


BAM_REGIONS = ["chr1", "chr2", "chr3", "chr1:1-1", "chr1:1-100", "chr1:500-1500", "chr1:9990-10110",
               "chr1:16384-16385", "chr1:15000-15000", "chr1:20000-20100", "chr1:16000-33000", "chr1:39900-40000",
               "chr2:1-50", "chr2:7050-7060", "chr2:12000-14000", "chr2:24950-25000", "chr2:24999-30000",
               "chr1:32768-32768", "chr1:5030-5050", "chr3:1-8000"]
DEPTH_REGIONS = ["chr1:1-3000", "chr1:9900-10200", "chr2:12000-14000", "chr2:24900-25000"]
TABIX_REGIONS = ["chr1", "chr2:1-1000", "chr1:1-1", "chr1:16380-16400", "chr1:5000-9000", "chr1:30000-40000",
                 "chr2:12345-12345", "chr2:20000-25000", "chr3", "chr3:7000-8000", "chr1:12000-12010"]


def write_regions(path, regions, cmd):
    with open(os.path.join(HERE, path), "w", newline="\n") as f:
        for reg in regions:
            f.write(f"## {reg}\n")
            text = capture(*cmd, reg).decode()
            for line in text.splitlines():
                c = line.split("\t")
                f.write("\t".join(c[:4]) + "\n")


def gzip_copy(src, dst):
    with open(os.path.join(HERE, src), "rb") as i, gzip.open(os.path.join(HERE, dst), "wb", 9) as o:
        o.write(i.read())


def generate():
    seqs = make_reference()
    make_sam(seqs)
    make_vcf(seqs)
    make_bed()
    make_gff()
    run("samtools", "faidx", "ref.fa")
    # BAM straight from the SAM text (record order kept) and the coordinate-sorted BAM + indexes
    run("samtools", "view", "--no-PG", "-b", "-o", "reads.unsorted.bam", "reads.sam")
    run("samtools", "sort", "--no-PG", "-o", "reads.bam", "reads.sam")
    run("samtools", "index", "reads.bam")
    run("samtools", "index", "-c", "reads.bam")
    with open(os.path.join(HERE, "reads.bam.sam"), "wb") as f:
        f.write(capture("samtools", "view", "--no-PG", "-h", "reads.bam"))
    gzip_copy("reads.bam.sam", "reads.bam.sam.gz")
    os.remove(os.path.join(HERE, "reads.bam.sam"))
    with open(os.path.join(HERE, "reads.idxstats.tsv"), "wb") as f:
        f.write(capture("samtools", "idxstats", "reads.bam"))
    write_regions("reads.regions.tsv", BAM_REGIONS, ["samtools", "view", "reads.bam"])
    # samtools depth with its defaults (exclude UNMAP,SECONDARY,QCFAIL,DUP; no deletions), and with -J -Q 30 -q 20
    with open(os.path.join(HERE, "reads.depth.tsv"), "w", newline="\n") as f:
        for args in [[], ["-J", "-Q", "30", "-q", "20"]]:
            for reg in DEPTH_REGIONS:
                f.write(f"## {' '.join(args)}|{reg}\n")
                f.write(capture("samtools", "depth", "-a", *args, "-r", reg, "reads.bam").decode())
    # CRAM 3.0 / 3.1 with different codec profiles, plus embedded / no reference and multi-reference slices
    crams = {
        "reads.30.cram": "version=3.0",
        "reads.31.cram": "version=3.1",
        "reads.31_small.cram": "version=3.1,small",
        "reads.31_archive.cram": "version=3.1,archive",
        "reads.30_bzip2.cram": "version=3.0,use_bzip2=1,use_rans=0",
        "reads.30_lzma.cram": "version=3.0,use_lzma=1,use_rans=0",
        "reads.30_embed.cram": "version=3.0,embed_ref=1",
        "reads.30_noref.cram": "version=3.0,no_ref=1",
        "reads.30_multi.cram": "version=3.0,seqs_per_slice=60,slices_per_container=3",
        "reads.31_lossy_names.cram": "version=3.1,lossy_names=1",
    }
    decoded, manifest = {}, []
    for name, opts in crams.items():
        run("samtools", "view", "--no-PG", "-C", *sum([["--output-fmt-option", o] for o in opts.split(",")], []), "-o", name, "reads.bam")
        run("samtools", "index", name)
        text = capture("samtools", "view", "--no-PG", "-h", name)
        # most profiles decode to identical SAM text: store each distinct expectation once
        expected = next((e for e, t in decoded.items() if t == text), None)
        if expected is None:
            expected = name + ".sam.gz"
            decoded[expected] = text
            with gzip.open(os.path.join(HERE, expected), "wb", 9) as f:
                f.write(text)
        manifest.append(f"{name}\t{opts}\t{expected}")
    with open(os.path.join(HERE, "cram.manifest.tsv"), "w", newline="\n") as f:
        f.write("\n".join(manifest) + "\n")
    write_regions("reads.30.cram.regions.tsv", BAM_REGIONS, ["samtools", "view", "reads.30.cram"])
    write_regions("reads.30_multi.cram.regions.tsv", BAM_REGIONS, ["samtools", "view", "reads.30_multi.cram"])
    # tabix: VCF (TBI + CSI), BED, GFF3
    for src, preset in [("variants.vcf", "vcf"), ("features.bed", "bed"), ("genes.gff3", "gff")]:
        with open(os.path.join(HERE, src + ".gz"), "wb") as f:
            f.write(capture("bgzip", "-c", src))
        run("tabix", "-f", "-p", preset, *(["-S", "1"] if preset == "bed" else []), src + ".gz")
        write_regions(src + ".regions.tsv", TABIX_REGIONS, ["tabix", src + ".gz"])
    run("tabix", "-f", "-C", "-p", "vcf", "variants.vcf.gz")
    gzip_copy("reads.sam", "reads.sam.gz")
    for f in ["variants.vcf", "features.bed", "genes.gff3", "reads.sam"]:
        os.remove(os.path.join(HERE, f))


# --------------------------------------------------------------------------------------------- hts-specs copy
SPEC_CRAM_SKIP = {"level-1.cram", "level-2.cram", "level-4.cram"}  # 0.5 MB each: too large for the repo


def copy_hts_specs(root):
    dst = os.path.join(HERE, "..", "hts-specs")
    src = os.path.join(root, "test", "cram")
    os.makedirs(os.path.join(dst, "cram30"), exist_ok=True)
    shutil.copy(os.path.join(src, "3.0", "LICENSE"), os.path.join(dst, "cram30", "LICENSE"))
    max_end = {}
    for name in sorted(os.listdir(os.path.join(src, "3.0", "passed"))):
        if name in SPEC_CRAM_SKIP:
            continue
        p = os.path.join(src, "3.0", "passed", name)
        if name.endswith(".sam"):
            with open(p) as f:
                text = f.read()
            for line in text.splitlines():
                if line.startswith("@") or not line:
                    continue
                c = line.split("\t")
                if c[2] != "*":
                    max_end[c[2]] = max(max_end.get(c[2], 0), int(c[3]) + max(len(c[9]), 1) + 2000)
            if len(text) > 20000:
                with gzip.open(os.path.join(dst, "cram30", name + ".gz"), "wt", newline="\n") as o:
                    o.write(text)
                continue
        shutil.copy(p, os.path.join(dst, "cram30", name))
    # ce.fa is 1 MB; keep only the leading bases the test reads touch. Slice MD5s cover only the slice span,
    # so they still verify against the truncated sequences.
    with open(os.path.join(src, "ce.fa")) as f:
        seqs, name = {}, None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    with open(os.path.join(dst, "cram30", "ce_trunc.fa"), "w", newline="\n") as o:
        for n, parts in seqs.items():
            s = "".join(parts)[:max(max_end.get(n, 0), 1000)]
            o.write(f">{n}\n")
            for i in range(0, len(s), 50):
                o.write(s[i:i + 50] + "\n")
    # codec vectors: a representative subset (all variants of the q4 / names / u32 data)
    cdst = os.path.join(dst, "codecs")
    for sub, keep in [("gzip", ["q4.gz", "q40+dir.gz", "01.names.gz", "20.names.gz", "u32.gz"]),
                      ("rans4x8", ["q4.0", "q4.1"]),
                      ("ransNx16", ["q4.0", "q4.1", "q4.4", "q4.5", "q4.64", "q4.65", "q4.128", "q4.129", "q4.192",
                                    "q4.193", "u32.1", "u32.9"]),
                      ("range", ["q4.0", "q4.1", "q4.8", "q4.9", "q4.64", "q4.65", "q4.128", "q4.129", "q4.192",
                                 "q4.193", "u32.1", "u32.9"]),
                      ("fqzcomp", ["q4.0", "q4.1", "q4.2", "q4.3", "q40+dir.3"]),
                      ("tok3", ["01.names.1", "01.names.3", "01.names.5", "01.names.7", "01.names.9", "01.names.11",
                                "01.names.13", "01.names.15", "01.names.17", "01.names.19", "20.names.1", "20.names.9",
                                "20.names.19"])]:
        os.makedirs(os.path.join(cdst, sub), exist_ok=True)
        for k in keep:
            shutil.copy(os.path.join(src, "codecs", sub, k), os.path.join(cdst, sub, k))
    shutil.copy(os.path.join(src, "codecs", "LICENSE"), os.path.join(cdst, "LICENSE"))
    shutil.copy(os.path.join(src, "codecs", "README.md"), os.path.join(cdst, "README.md"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hts-specs", help="hts-specs checkout to copy conformance files from")
    ap.add_argument("--skip-generate", action="store_true")
    a = ap.parse_args()
    if not a.skip_generate:
        generate()
    if a.hts_specs:
        copy_hts_specs(a.hts_specs)
