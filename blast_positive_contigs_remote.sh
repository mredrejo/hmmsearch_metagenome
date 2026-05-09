#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
USO:
  ./blast_positive_contigs_remote.sh -i results_pipeline -o blast_species [opciones]

DESCRIPCIÓN:
  A partir de la carpeta de resultados del pipeline MAG/HMM, extrae los contigs positivos
  presentes en hits.tsv, ejecuta blastn --remote para cada contig y genera:
    1) Un FASTA con todos los contigs positivos.
    2) Un BLAST por contig en formato outfmt 7.
    3) Una tabla resumen con el top hit de cada contig y el SRA/ID de procedencia.

REQUISITOS DE ENTRADA:
  results_pipeline/hits.tsv
  results_pipeline/nt/*.fna, normalmente:
    - SRA: nt/SRRxxxx.contigs.fna, nt/ERRxxxx.contigs.fna o nt/DRRxxxx.contigs.fna
    - WGS: nt/<ID>.fna o similar

OPCIONES:
  -i, --indir DIR          Carpeta de resultados del pipeline original (requerido)
  -o, --outdir DIR         Carpeta de salida para BLAST (default: blast_species)
  --db DB                  Base de datos BLAST remota (default: nt)
  --entrez-query QUERY     Filtro Entrez opcional, ej. 'bacteria[organism]'
  --max-target-seqs N      Máximo de hits por contig (default: 10)
  --evalue X               E-value para blastn (default: 1e-20)
  --word-size N            Word size para blastn (default: 28)
  --perc-identity X        Identidad mínima opcional, ej. 80. Si se omite no se usa
  --sleep SEC              Pausa entre llamadas remotas a NCBI (default: 5)
  --resume                 No repetir BLAST si ya existe un .blast7 no vacío
  -h, --help               Muestra esta ayuda

SALIDAS:
  outdir/positive_contigs.fna
  outdir/contig_map.tsv
  outdir/blast7/<contig>.blast7
  outdir/blast_top_hits.tsv
  outdir/logs/blast_species.log

NOTAS:
  - Este script usa blastn --remote, por lo que requiere conexión a internet.
  - Para evitar nombres problemáticos en ficheros, los IDs de contig se sanean en los .blast7.
  - El resumen usa el primer hit no comentado de cada archivo outfmt 7 como top hit.
EOF
}

INDIR=""
OUTDIR="blast_species"
DB="nt"
ENTREZ_QUERY=""
MAX_TARGET_SEQS=10
EVALUE="1e-20"
WORD_SIZE=28
PERC_IDENTITY=""
SLEEP_SEC=5
RESUME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir) INDIR="${2:-}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:-}"; shift 2 ;;
    --db) DB="${2:-}"; shift 2 ;;
    --entrez-query) ENTREZ_QUERY="${2:-}"; shift 2 ;;
    --max-target-seqs) MAX_TARGET_SEQS="${2:-}"; shift 2 ;;
    --evalue) EVALUE="${2:-}"; shift 2 ;;
    --word-size) WORD_SIZE="${2:-}"; shift 2 ;;
    --perc-identity) PERC_IDENTITY="${2:-}"; shift 2 ;;
    --sleep) SLEEP_SEC="${2:-}"; shift 2 ;;
    --resume) RESUME=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Opción desconocida: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INDIR" ]] || { echo "[ERROR] Falta -i/--indir" >&2; usage; exit 1; }
[[ -d "$INDIR" ]] || { echo "[ERROR] No existe carpeta: $INDIR" >&2; exit 1; }
[[ -s "$INDIR/hits.tsv" ]] || { echo "[ERROR] No existe o está vacío: $INDIR/hits.tsv" >&2; exit 1; }
[[ -d "$INDIR/nt" ]] || { echo "[ERROR] No existe carpeta de contigs: $INDIR/nt" >&2; exit 1; }

for cmd in blastn python3 awk grep sed sort; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] No encuentro dependencia: $cmd" >&2; exit 1; }
done

mkdir -p "$OUTDIR" "$OUTDIR/blast7" "$OUTDIR/tmp" "$OUTDIR/logs"
LOG="$OUTDIR/logs/blast_species.log"
: > "$LOG"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2
}

POS_FASTA="$OUTDIR/positive_contigs.fna"
MAP_TSV="$OUTDIR/contig_map.tsv"
SUMMARY="$OUTDIR/blast_top_hits.tsv"

log "Inicio. INDIR=$INDIR OUTDIR=$OUTDIR DB=$DB EVALUE=$EVALUE MAX_TARGET_SEQS=$MAX_TARGET_SEQS"

python3 - "$INDIR" "$POS_FASTA" "$MAP_TSV" <<'PY'
import sys, csv, re
from pathlib import Path

indir = Path(sys.argv[1])
out_fasta = Path(sys.argv[2])
map_tsv = Path(sys.argv[3])
hits = indir / "hits.tsv"
nt_dir = indir / "nt"

# Leer hits.tsv y detectar columnas de forma tolerante.
with hits.open(newline="") as fh:
    reader = csv.reader(fh, delimiter="\t")
    rows = list(reader)

if not rows:
    raise SystemExit("hits.tsv vacío")

header = rows[0]
data = rows[1:] if any(x.lower() in ("contig_id", "sra", "sample", "run") for x in header) else rows
header_lower = [h.lower() for h in header]

def find_col(candidates, default=None):
    for c in candidates:
        if c in header_lower:
            return header_lower.index(c)
    return default

contig_col = find_col(["contig_id", "contig", "target", "target_name"], 0)
sra_col = find_col(["sra", "run", "sample", "sample_id", "source", "assembly", "master"], None)

positive = []
seen = set()
for row in data:
    if not row or len(row) <= contig_col:
        continue
    contig = row[contig_col].strip()
    if not contig or contig.lower() == "contig_id":
        continue
    # Inferir SRA/ID desde columna si existe; si no, desde contig o más adelante desde el fichero.
    sra = ""
    if sra_col is not None and len(row) > sra_col:
        sra = row[sra_col].strip()
    if not sra:
        m = re.search(r"\b([SED]RR\d+)\b", contig)
        sra = m.group(1) if m else "NA"
    key = (contig, sra)
    if key not in seen:
        positive.append((contig, sra))
        seen.add(key)

if not positive:
    raise SystemExit("No se detectaron contigs positivos en hits.tsv")

# Indexar FASTA disponibles en nt/. Para ahorrar memoria, se hace un índice simple contig -> fichero.
fasta_files = sorted(list(nt_dir.glob("*.fna")) + list(nt_dir.glob("*.fa")) + list(nt_dir.glob("*.fasta")))
if not fasta_files:
    raise SystemExit(f"No hay FASTA en {nt_dir}")

wanted = {c for c, s in positive}
found = {}
seq_chunks = {}
current_id = None
current_header = None
current_chunks = []

def flush_record(file_path):
    global current_id, current_header, current_chunks
    if current_id in wanted and current_id not in found:
        found[current_id] = (file_path, current_header, "".join(current_chunks))

for fp in fasta_files:
    current_id = None
    current_header = None
    current_chunks = []
    with fp.open(errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    flush_record(fp)
                current_header = line[1:]
                current_id = current_header.split()[0]
                current_chunks = []
            else:
                if current_id is not None:
                    current_chunks.append(line.strip())
        if current_id is not None:
            flush_record(fp)

missing = [(c,s) for c,s in positive if c not in found]

def sanitize_id(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", x)

with out_fasta.open("w") as out, map_tsv.open("w") as mp:
    mp.write("query_id\toriginal_contig_id\tsource_id\tcontig_file\tlength\n")
    for contig, sra in positive:
        if contig not in found:
            continue
        fp, header, seq = found[contig]
        # Si no hay SRA/ID en hits, intentar inferir del nombre del FASTA: SRR.contigs.fna -> SRR
        source = sra
        if source == "NA" or not source:
            stem = fp.name
            m = re.search(r"\b([SED]RR\d+)\b", stem)
            source = m.group(1) if m else stem.split(".")[0]
        qid = sanitize_id(contig)
        out.write(f">{qid} original_contig_id={contig} source_id={source} contig_file={fp.name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")
        mp.write(f"{qid}\t{contig}\t{source}\t{fp.name}\t{len(seq)}\n")

if missing:
    miss_file = out_fasta.parent / "missing_positive_contigs.tsv"
    with miss_file.open("w") as mfh:
        mfh.write("contig_id\tsource_id\n")
        for c,s in missing:
            mfh.write(f"{c}\t{s}\n")
    print(f"WARNING: {len(missing)} contigs positivos no se encontraron en nt/*.fna. Ver {miss_file}", file=sys.stderr)

print(f"Extracted {len(positive)-len(missing)} positive contigs to {out_fasta}", file=sys.stderr)
PY

n_contigs=$(grep -c '^>' "$POS_FASTA" || true)
log "Contigs positivos extraídos: $n_contigs"
[[ "$n_contigs" -gt 0 ]] || { log "ERROR: no se extrajo ningún contig positivo"; exit 1; }

# Ejecutar blastn remoto contig por contig.
# outfmt 7 con campos explícitos para poder parsear una tabla resumen robusta.
OUTFMT='7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sscinames scomnames staxids'

awk '/^>/{if(seq){print id"\t"seq}; id=substr($1,2); seq=""; next} {seq=seq$0} END{if(seq){print id"\t"seq}}' "$POS_FASTA" \
| while IFS=$'\t' read -r qid seq; do
    qfa="$OUTDIR/tmp/${qid}.fa"
    bout="$OUTDIR/blast7/${qid}.blast7"
    printf '>%s\n%s\n' "$qid" "$seq" > "$qfa"

    if [[ "$RESUME" -eq 1 && -s "$bout" ]]; then
      log "SKIP $qid: ya existe $bout"
      continue
    fi

    log "BLAST remoto para $qid"
    args=(blastn -query "$qfa" -db "$DB" -remote -outfmt "$OUTFMT" -max_target_seqs "$MAX_TARGET_SEQS" -evalue "$EVALUE" -word_size "$WORD_SIZE" -out "$bout")
    [[ -n "$ENTREZ_QUERY" ]] && args+=(-entrez_query "$ENTREZ_QUERY")
    [[ -n "$PERC_IDENTITY" ]] && args+=(-perc_identity "$PERC_IDENTITY")

    if ! "${args[@]}" >> "$LOG" 2>&1; then
      log "ERROR: blastn falló para $qid. Revisa $LOG"
      # Crear archivo marcador para que quede constancia del fallo.
      {
        echo "# BLASTN failed for $qid"
        echo "# See $LOG"
      } > "$bout"
    fi

    sleep "$SLEEP_SEC"
  done

# Resumen: primer hit no comentado de cada .blast7.
python3 - "$MAP_TSV" "$OUTDIR/blast7" "$SUMMARY" <<'PY'
import sys, csv
from pathlib import Path

map_tsv = Path(sys.argv[1])
blast_dir = Path(sys.argv[2])
summary = Path(sys.argv[3])

meta = {}
with map_tsv.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        meta[row["query_id"]] = row

fields = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle","sscinames","scomnames","staxids"]

with summary.open("w") as out:
    out.write("query_id\toriginal_contig_id\tsource_id\tcontig_file\tcontig_length\tstatus\t" + "\t".join(fields[1:]) + "\n")
    for qid in sorted(meta):
        row = meta[qid]
        blast_file = blast_dir / f"{qid}.blast7"
        top = None
        status = "NO_BLAST_FILE"
        if blast_file.exists():
            status = "NO_HIT"
            with blast_file.open(errors="replace") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= len(fields):
                        top = dict(zip(fields, parts[:len(fields)]))
                        status = "OK"
                        break
                    else:
                        status = "PARSE_ERROR"
                        break
        base = [qid, row.get("original_contig_id",""), row.get("source_id",""), row.get("contig_file",""), row.get("length",""), status]
        if top:
            vals = [top.get(f, "") for f in fields[1:]]
        else:
            vals = ["NA"] * (len(fields)-1)
        out.write("\t".join(base + vals) + "\n")
PY

log "Resumen final creado: $SUMMARY"
log "Fin"
