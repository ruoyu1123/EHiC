#!/usr/bin/env python3
import gzip
import json
import math
import re
import struct
import sys
import zlib
from pathlib import Path


def read_png(path):
    data = Path(path).read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"not a PNG: {path}")

    pos = 8
    width = height = color_type = bit_depth = None
    chunks = []
    while pos < len(data):
        length = struct.unpack(">I", data[pos:pos + 4])[0]
        kind = data[pos + 4:pos + 8]
        payload = data[pos + 8:pos + 8 + length]
        pos += 12 + length
        if kind == b"IHDR":
            width, height, bit_depth, color_type, _, _, _ = struct.unpack(">IIBBBBB", payload)
            if bit_depth != 8 or color_type not in (0, 2, 6):
                raise ValueError("only 8-bit grayscale/RGB/RGBA PNGs are supported")
        elif kind == b"IDAT":
            chunks.append(payload)
        elif kind == b"IEND":
            break

    channels = {0: 1, 2: 3, 6: 4}[color_type]
    raw = zlib.decompress(b"".join(chunks))
    stride = width * channels
    rows = []
    prev = [0] * stride
    i = 0
    for _ in range(height):
        filt = raw[i]
        i += 1
        row = list(raw[i:i + stride])
        i += stride
        recon = [0] * stride
        for x, value in enumerate(row):
            left = recon[x - channels] if x >= channels else 0
            up = prev[x]
            up_left = prev[x - channels] if x >= channels else 0
            if filt == 0:
                pred = 0
            elif filt == 1:
                pred = left
            elif filt == 2:
                pred = up
            elif filt == 3:
                pred = (left + up) // 2
            elif filt == 4:
                p = left + up - up_left
                pa = abs(p - left)
                pb = abs(p - up)
                pc = abs(p - up_left)
                pred = left if pa <= pb and pa <= pc else up if pb <= pc else up_left
            else:
                raise ValueError(f"unknown PNG filter {filt}")
            recon[x] = (value + pred) & 255
        prev = recon
        pixels = []
        for x in range(width):
            base = x * channels
            if color_type == 0:
                gray = recon[base]
                pixels.append((gray, gray, gray))
            else:
                pixels.append(tuple(recon[base:base + 3]))
        rows.append(pixels)
    return width, height, rows


def write_png_rgb(path, pixels):
    height = len(pixels)
    width = len(pixels[0]) if height else 0
    raw = bytearray()
    for row in pixels:
        raw.append(0)
        for r, g, b in row:
            raw.extend((r & 255, g & 255, b & 255))
    def chunk(kind, payload):
        return (struct.pack(">I", len(payload)) + kind + payload +
                struct.pack(">I", zlib.crc32(kind + payload) & 0xffffffff))
    png = bytearray(b"\x89PNG\r\n\x1a\n")
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(bytes(raw), 9))
    png += chunk(b"IEND", b"")
    Path(path).write_bytes(png)


def resize_average(matrix, size):
    src_h = len(matrix)
    src_w = len(matrix[0])
    out = [[0.0] * size for _ in range(size)]
    for y in range(size):
        y0 = int(y * src_h / size)
        y1 = max(y0 + 1, int((y + 1) * src_h / size))
        for x in range(size):
            x0 = int(x * src_w / size)
            x1 = max(x0 + 1, int((x + 1) * src_w / size))
            total = count = 0.0
            for yy in range(y0, min(y1, src_h)):
                for xx in range(x0, min(x1, src_w)):
                    total += matrix[yy][xx]
                    count += 1.0
            out[y][x] = total / count if count else 0.0
    return out


def extract_contact_image(path, bins):
    width, height, pixels = read_png(path)
    colored = []
    for y, row in enumerate(pixels):
        for x, (r, g, b) in enumerate(row):
            saturation = max(r, g, b) - min(r, g, b)
            redish = r > g + 8 and r > b + 8
            dark = r < 90 and g < 90 and b < 90
            if saturation > 18 or redish or dark:
                colored.append((x, y))
    if not colored:
        raise ValueError("could not locate heatmap body in bacteria.png")
    min_x = min(x for x, _ in colored)
    max_x = max(x for x, _ in colored)
    min_y = min(y for _, y in colored)
    max_y = max(y for _, y in colored)
    side = min(max_x - min_x + 1, max_y - min_y + 1)
    cx = (min_x + max_x) // 2
    cy = (min_y + max_y) // 2
    left = max(0, min(width - side, cx - side // 2))
    top = max(0, min(height - side, cy - side // 2))

    body = []
    for y in range(top, top + side):
        line = []
        for x in range(left, left + side):
            r, g, b = pixels[y][x]
            luma = 0.299 * r + 0.587 * g + 0.114 * b
            red_signal = max(0.0, (r - 0.5 * (g + b)) / 255.0)
            darkness = max(0.0, 1.0 - luma / 255.0)
            line.append(max(red_signal, darkness))
        body.append(line)

    resized = resize_average(body, bins)
    for i in range(bins):
        for j in range(i, bins):
            v = 0.5 * (resized[i][j] + resized[j][i])
            resized[i][j] = resized[j][i] = v
    return resized, {"png_width": width, "png_height": height, "crop_left": left, "crop_top": top, "crop_side": side}


def write_sparse_matrix(path, matrix, min_weight=1e-6):
    values = [v for row in matrix for v in row if v > 0.0]
    floor = sorted(values)[int(0.05 * (len(values) - 1))] if values else 0.0
    scale = 1000.0
    with Path(path).open("w", newline="\n") as out:
        out.write("bin1\tbin2\tvalue\n")
        for i, row in enumerate(matrix):
            for j in range(i, len(row)):
                v = max(0.0, row[j] - floor)
                distance = abs(i - j)
                if distance == 0:
                    v *= 2.5
                elif distance <= 3:
                    v *= 1.25
                weight = v * scale
                if weight > min_weight:
                    out.write(f"{i}\t{j}\t{weight:.6f}\n")


def decompress_fasta(gz_path, out_path):
    with gzip.open(gz_path, "rt") as src, Path(out_path).open("w", newline="\n") as dst:
        for line in src:
            dst.write(line)


def fasta_name_and_length(path):
    name = None
    total = 0
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].split()[0]
                continue
            if name is not None:
                total += len(line)
    if not name or total == 0:
        raise ValueError(f"empty FASTA: {path}")
    return name, total


HEADER_RE = re.compile(r"^@hicreate_\d+_bin(\d+)_bin(\d+)/1$")


def matrix_from_fastq_headers(path, bins):
    counts = [[0.0] * bins for _ in range(bins)]
    pairs = 0
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 != 0:
                continue
            match = HEADER_RE.match(line.strip())
            if not match:
                continue
            i = int(match.group(1))
            j = int(match.group(2))
            if i < bins and j < bins:
                counts[i][j] += 1.0
                if i != j:
                    counts[j][i] += 1.0
                pairs += 1
    return counts, pairs


def scale_matrix_to_size(matrix, size):
    return resize_average(matrix, size)


def log_normalize(matrix):
    logged = [[math.log1p(v) for v in row] for row in matrix]
    max_v = max(max(row) for row in logged) if logged else 0.0
    if max_v <= 0.0:
        return logged
    return [[v / max_v for v in row] for row in logged]


def pearson(a, b):
    va = [x for row in a for x in row]
    vb = [x for row in b for x in row]
    n = min(len(va), len(vb))
    va = va[:n]
    vb = vb[:n]
    ma = sum(va) / n
    mb = sum(vb) / n
    num = sum((va[i] - ma) * (vb[i] - mb) for i in range(n))
    da = math.sqrt(sum((x - ma) ** 2 for x in va))
    db = math.sqrt(sum((x - mb) ** 2 for x in vb))
    return num / (da * db) if da > 0.0 and db > 0.0 else 0.0


def mse(a, b):
    values = [(a[y][x] - b[y][x]) ** 2 for y in range(len(a)) for x in range(len(a[y]))]
    return sum(values) / len(values) if values else 0.0


def red_colormap(v):
    v = max(0.0, min(1.0, v))
    if v < 0.7:
        t = v / 0.7
        r = 255
        g = int(255 * (1.0 - 0.78 * t))
        b = int(255 * (1.0 - 0.90 * t))
    else:
        t = (v - 0.7) / 0.3
        r = int(255 * (1.0 - t))
        g = int(56 * (1.0 - t))
        b = int(26 * (1.0 - t))
    return r, g, b


def blue_red_colormap(v):
    v = max(-1.0, min(1.0, v))
    if v < 0:
        t = -v
        return int(255 * (1.0 - t)), int(255 * (1.0 - t)), 255
    t = v
    return 255, int(255 * (1.0 - t)), int(255 * (1.0 - t))


def render_heatmap(matrix, path, cell=4):
    norm = log_normalize(matrix)
    pixels = []
    for row in norm:
        row_pixels = []
        for v in row:
            row_pixels.extend([red_colormap(v)] * cell)
        for _ in range(cell):
            pixels.append(list(row_pixels))
    write_png_rgb(path, pixels)


def render_difference(a, b, path, cell=4):
    pixels = []
    for y in range(len(a)):
        row_pixels = []
        for x in range(len(a[y])):
            row_pixels.extend([blue_red_colormap(a[y][x] - b[y][x])] * cell)
        for _ in range(cell):
            pixels.append(list(row_pixels))
    write_png_rgb(path, pixels)


def main(argv):
    if len(argv) < 2:
        print("usage: bacteria_hic_pipeline.py prepare|compare ...", file=sys.stderr)
        return 2
    mode = argv[1]
    if mode == "prepare":
        image = argv[2]
        gz_fasta = argv[3]
        out_prefix = Path(argv[4])
        bins = int(argv[5])
        matrix, crop = extract_contact_image(image, bins)
        matrix_path = f"{out_prefix}_matrix.tsv"
        source_heatmap = f"{out_prefix}_source_4m_heatmap.png"
        write_sparse_matrix(matrix_path, matrix)
        render_heatmap(matrix, source_heatmap)
        fasta_path = f"{out_prefix}_reference.fna"
        decompress_fasta(gz_fasta, fasta_path)
        contig, length = fasta_name_and_length(fasta_path)
        offset_path = f"{out_prefix}_offset.tsv"
        with Path(offset_path).open("w", newline="\n") as out:
            out.write("contig\tstart_bin\tend_bin\n")
            out.write(f"{contig}\t0\t{bins}\n")
        meta = {
            "source_image": image,
            "reference_fasta": fasta_path,
            "reference_contig": contig,
            "reference_length_bp": length,
            "simulated_source_length_bp": 4_000_000,
            "bin_size_bp": 20_000,
            "source_bins": bins,
            "matrix": matrix_path,
            "offset": offset_path,
            "source_heatmap": source_heatmap,
            "crop": crop,
        }
        Path(f"{out_prefix}_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        print(json.dumps(meta, indent=2))
        return 0
    if mode == "compare":
        image = argv[2]
        fastq_r1 = argv[3]
        out_prefix = Path(argv[4])
        source_bins = int(argv[5])
        target_bins = int(argv[6])
        source, crop = extract_contact_image(image, source_bins)
        source_target = log_normalize(scale_matrix_to_size(source, target_bins))
        reads, pairs = matrix_from_fastq_headers(fastq_r1, target_bins)
        reads_norm = log_normalize(reads)
        render_heatmap(reads, f"{out_prefix}_reads_heatmap.png")
        render_heatmap(scale_matrix_to_size(source, target_bins), f"{out_prefix}_source_rescaled_heatmap.png")
        render_difference(reads_norm, source_target, f"{out_prefix}_difference_reads_minus_source.png")
        stats = {
            "pairs_used": pairs,
            "source_bins": source_bins,
            "target_bins": target_bins,
            "pearson_log_normalized": pearson(reads_norm, source_target),
            "mse_log_normalized": mse(reads_norm, source_target),
            "reads_heatmap": f"{out_prefix}_reads_heatmap.png",
            "source_rescaled_heatmap": f"{out_prefix}_source_rescaled_heatmap.png",
            "difference_heatmap": f"{out_prefix}_difference_reads_minus_source.png",
            "crop": crop,
        }
        Path(f"{out_prefix}_compare.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
        print(json.dumps(stats, indent=2))
        return 0
    print(f"unknown mode: {mode}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
