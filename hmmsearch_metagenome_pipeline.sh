#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
USO:

Modo WGS:
  ./mag_pipeline_v2_unified_nocompress_noflock.sh --wgs -i masters.tsv -p profile.hmm [opciones]

Modo SRA desde BioProject:
  ./mag_pipeline_v2_unified_nocompress_noflock.sh --bioproject PRJNAxxxxxx -p profile.hmm [opciones]

Modo SRA desde lista de SRR/ERR/DRR:
  ./mag_pipeline_v2_unified_nocompress_noflock.sh --sra-list sra_ids.txt -p profile.hmm [opciones]

Opciones comunes:
  -p, --profile FILE        Perfil HMM (requerido)
  -o, --outdir DIR          Directorio salida (default: results_pipeline)
  --evalue X                Umbral E-value por SECUENCIA (default: 1e-5)
  --cpu N                   CPUs hmmsearch (default: 4)

Paralelismo:
  --jobs N                  Jobs por defecto (default: 2)
  --jobs-wgs N              Jobs WGS (default: jobs)
  --jobs-sra N              Jobs SRA (default: jobs)

WGS:
  --wgs
  -i, --input FILE          TSV/CSV con masters
  --col N                   Columna con ID master (default: 1)
  --keep-negatives          No borrar nt/aa/hmm de negativos (default: borra)
  --wgs-max-parts N         Máximo de partes wgs.PROJ.1..N (default: 200)

SRA:
  --bioproject PRJNA...      Obtiene SRR automáticamente desde BioProject
  --sra-list FILE            TXT con una lista de SRR/ERR/DRR, uno por línea
  --max-runs N              Limita SRR (default: 0=todos)
  --sra-threads N           Threads fasterq-dump (default: 6)
  --assembler metaspades|spades_meta (default: metaspades)
  --spades-threads N        Threads por SRR (default: 8)
  --spades-memory GB        RAM por SRR (default: 32)
  --keep-assembly           Conserva assembly aunque negativo
  --keep-reads              Conserva reads (default: borra)

Salidas:
  outdir/hits.tsv
  outdir/summary.tsv
  outdir/logs/*.log
  outdir/logs/pipeline.log          Log general de ejecución
  outdir/logs/<SRR>.spades.log      Log nativo de SPAdes si existe y SPAdes falla
EOF
}

# ---------------- defaults ----------------
MODE=""
INPUT=""
COL=1
BIOPROJECT=""
SRA_LIST=""
PROFILE=""
OUTDIR="results_pipeline"
EVALUE="1e-5"
CPU=4

JOBS=2
JOBS_WGS=""
JOBS_SRA=""

KEEP_NEG=0
WGS_MAX_PARTS=200

MAX_RUNS=0
SRA_THREADS=6
ASSEMBLER="metaspades"
SPADES_THREADS=8
SPADES_MEM=32
KEEP_ASM=0
KEEP_READS=0

TMPDIR=""

# ---------------- parse args ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --wgs) MODE="WGS"; shift ;;
    --bioproject) MODE="SRA"; BIOPROJECT="${2:-}"; shift 2 ;;
    --sra-list) MODE="SRA"; SRA_LIST="${2:-}"; shift 2 ;;
    -i|--input) INPUT="${2:-}"; shift 2 ;;
    --col) COL="${2:-}"; shift 2 ;;
    -p|--profile) PROFILE="${2:-}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:-}"; shift 2 ;;
    --evalue) EVALUE="${2:-}"; shift 2 ;;
    --cpu) CPU="${2:-}"; shift 2 ;;
    --jobs) JOBS="${2:-}"; shift 2 ;;
    --jobs-wgs) JOBS_WGS="${2:-}"; shift 2 ;;
    --jobs-sra) JOBS_SRA="${2:-}"; shift 2 ;;
    --keep-negatives) KEEP_NEG=1; shift ;;
    --wgs-max-parts) WGS_MAX_PARTS="${2:-}"; shift 2 ;;
    --max-runs) MAX_RUNS="${2:-}"; shift 2 ;;
    --sra-threads) SRA_THREADS="${2:-}"; shift 2 ;;
    --assembler) ASSEMBLER="${2:-}"; shift 2 ;;
    --spades-threads) SPADES_THREADS="${2:-}"; shift 2 ;;
    --spades-memory) SPADES_MEM="${2:-}"; shift 2 ;;
    --keep-assembly) KEEP_ASM=1; shift ;;
    --keep-reads) KEEP_READS=1; shift ;;
    --tmpdir) TMPDIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Opción desconocida: $1"; usage; exit 1 ;;
  esac
done

[[ -n "$MODE" ]] || { echo "[ERROR] Debes indicar --wgs, --bioproject o --sra-list"; usage; exit 1; }
[[ -f "$PROFILE" ]] || { echo "[ERROR] No existe perfil: $PROFILE"; exit 1; }

JOBS_WGS="${JOBS_WGS:-$JOBS}"
JOBS_SRA="${JOBS_SRA:-$JOBS}"

mkdir -p "$OUTDIR"/{nt,aa,hmm,logs,tmp,assembly,reads,parts}
TMPDIR="${TMPDIR:-$OUTDIR/tmp}"
mkdir -p "$TMPDIR"

PIPELINE_LOG="$OUTDIR/logs/pipeline.log"
: > "$PIPELINE_LOG"

log_msg() {
  # Log general. Incluye fecha ISO-8601 y PID para distinguir jobs en paralelo.
  # No se usa flock para mantener esta versión compatible con el nombre noflock.
  local level="$1"; shift
  printf '[%s] [%s] [pid:%s] %s
' "$(date '+%Y-%m-%d %H:%M:%S')" "$level" "$$" "$*" >> "$PIPELINE_LOG"
}

log_info() { log_msg "INFO" "$@"; }
log_warn() { log_msg "WARN" "$@"; }
log_error() { log_msg "ERROR" "$@"; }

HITS="$OUTDIR/hits.tsv"
SUMMARY="$OUTDIR/summary.tsv"

HITS_HDR="contig_id\tprodigal_protein_id\thmm_query\tseq_evalue\tseq_score\tseq_bias\tsource_id"
SUM_HDR="source_id\tstatus\thits\tbest_seq_evalue\tbest_seq_score"

# ---------------- deps ----------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Falta comando: $1"; exit 1; }; }
for cmd in awk sed grep sort prodigal hmmsearch curl gunzip efetch; do need "$cmd"; done
if [[ "$MODE" == "SRA" ]]; then
  for cmd in prefetch fasterq-dump; do need "$cmd"; done
  # esearch/efetch solo son necesarios cuando los SRR se obtienen desde BioProject.
  if [[ -n "$BIOPROJECT" ]]; then
    for cmd in esearch efetch; do need "$cmd"; done
  fi
  if [[ "$ASSEMBLER" == "metaspades" ]]; then
    need metaspades.py
  else
    need spades.py
  fi
fi

# ---------------- helpers ----------------
normalize_fasta_headers() {
  local infa="$1" outfa="$2"
  awk '
    /^>/{
      sub(/^>/,"",$0);
      split($0,a,/[\t ]+/);
      print ">" a[1];
      next
    }
    {print}
  ' "$infa" > "$outfa"
}

summarize_tbl() {
  local tbl="$1"
  awk '
    BEGIN{FS="[ \t]+"; hits=0; bestE=""; bestS=""}
    /^#/ {next}
    NF>0 {
      hits++; e=$5; s=$6;
      if(bestE=="" || (e+0 < bestE+0)) { bestE=e; bestS=s }
    }
    END{
      if(hits==0) print "0\tNA\tNA";
      else print hits "\t" bestE "\t" bestS
    }
  ' "$tbl"
}

tblout_to_hits_sra() {
  # source_id constante = SRR
  local tbl="$1" srr="$2"
  awk -v SRC="$srr" '
    BEGIN{FS="[ \t]+"; OFS="\t"}
    /^#/ {next}
    NF>0 {
      prot=$1; hmm=$3; e=$5; s=$6; b=$7;
      contig=prot; sub(/_[0-9]+$/, "", contig);
      print contig, prot, hmm, e, s, b, SRC
    }
  ' "$tbl"
}

tblout_to_hits_wgs() {
  # source_id = contig_id (NO master, NO token)
  local tbl="$1"
  awk '
    BEGIN{FS="[ \t]+"; OFS="\t"}
    /^#/ {next}
    NF>0 {
      prot=$1; hmm=$3; e=$5; s=$6; b=$7;
      contig=prot; sub(/_[0-9]+$/, "", contig);
      print contig, prot, hmm, e, s, b, contig
    }
  ' "$tbl"
}

pick_read_file() {
  local base="$1"
  if [[ -f "${base}.fastq.gz" ]]; then echo "${base}.fastq.gz"
  elif [[ -f "${base}.fastq" ]]; then echo "${base}.fastq"
  else echo ""
  fi
}

# ---------------- simple parallel queue ----------------
run_queue() {
  local jobs_max="$1"; shift
  local listfile="$1"; shift
  local worker_fn="$1"; shift

  local -a pids=()

  while IFS= read -r item; do
    [[ -z "$item" ]] && continue

    # lanzar job
    "$worker_fn" "$item" &
    pids+=( "$!" )

    # limitar concurrencia: cuando llegamos al máximo, esperamos al más antiguo
    if (( ${#pids[@]} >= jobs_max )); then
      wait "${pids[0]}" || true
      pids=( "${pids[@]:1}" )
    fi
  done < "$listfile"

  # esperar los restantes
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || true
  done
}

# ---------------- WGS workers ----------------
extract_master_id() {
  local line="$1"
  local id
  id="$(echo "$line" | awk -v c="$COL" 'BEGIN{FS="\t"} {if(NF>=c){print $c; exit}}')"
  if [[ -z "$id" ]]; then
    id="$(echo "$line" | awk -v c="$COL" 'BEGIN{FS=","} {if(NF>=c){print $c; exit}}')"
  fi
  echo "${id//[[:space:]]/}"
}

get_wgs_token_from_master() {
  local master="$1"
  efetch -db nuccore -id "$master" -format gb | awk '
    /^[[:space:]]*(WGS|TSA)[[:space:]]+/ {
      line=$0
      sub(/^[[:space:]]*/, "", line)
      sub(/^(WGS|TSA)[[:space:]]+/, "", line)

      # Primer token tras WGS/TSA (p.ej. JNMI01000001-JNMI01000178)
      split(line, a, /[ \t]+/)
      x=a[1]

      # Si es rango, quedarnos con la parte izquierda (JNMI01000001)
      split(x, b, /-/)
      left=b[1]

      # Extraer token: 4-6 letras + 2 dígitos, SIN usar {m,n}
      # ^[A-Z][A-Z][A-Z][A-Z][A-Z]?[A-Z]?[0-9][0-9]
      if (match(left, /^[A-Z][A-Z][A-Z][A-Z][A-Z]?[A-Z]?[0-9][0-9]/)) {
        print substr(left, RSTART, RLENGTH)
        exit
      }
    }
  '
}


wgs_subdir_for_project() {
  local proj="$1"
  if [[ ${#proj} -eq 4 ]]; then echo "${proj:0:1}"; else echo "${proj:0:3}"; fi
}

download_wgs_fasta_nt_ftp() {
  local proj="$1" outfa="$2" log="$3" maxparts="$4"
  local subdir; subdir="$(wgs_subdir_for_project "$proj")"
  : > "$outfa"
  local ok_any=0

  for part in $(seq 1 "$maxparts"); do
    local url="https://ftp.ncbi.nlm.nih.gov/genbank/wgs/${subdir}/wgs.${proj}.${part}.fsa_nt.gz"
    echo "[INFO] FTP: $url" >> "$log"
    if curl -fSL "$url" 2>>"$log" | gunzip -c >> "$outfa" 2>>"$log"; then
      ok_any=1
    else
      break
    fi
  done

  [[ "$ok_any" -eq 1 ]] && [[ -s "$outfa" ]] && head -n 1 "$outfa" | grep -q '^>' && return 0
  return 1
}

process_one_wgs_master() {
  local master="$1"
  local log="$OUTDIR/logs/${master}.wgs.log"
  : > "$log"

  local token; token="$(get_wgs_token_from_master "$master")"
  if [[ -z "$token" ]]; then
    echo -e "${master}\tWGS_TOKEN_FAIL\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    return 0
  fi

  local proj; proj="$(echo "$token" | sed -nE 's/^([A-Z]{4,6}).*/\1/p')"
  if [[ -z "$proj" ]]; then
    echo -e "${master}\tWGS_CODE_FAIL\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    return 0
  fi

  local raw="$OUTDIR/nt/${master}.raw.fna"
  local ntfa="$OUTDIR/nt/${master}.fna"
  local aafa="$OUTDIR/aa/${master}.faa"
  local gbk="$OUTDIR/nt/${master}.prodigal.gbk"
  local tbl="$OUTDIR/hmm/${master}.tblout"

  if ! download_wgs_fasta_nt_ftp "$proj" "$raw" "$log" "$WGS_MAX_PARTS"; then
    echo -e "${master}\tWGS_DOWNLOAD_FAIL\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    rm -f "$raw" 2>/dev/null || true
    return 0
  fi

  normalize_fasta_headers "$raw" "$ntfa"
  rm -f "$raw"

  if ! prodigal -i "$ntfa" -a "$aafa" -o "$gbk" -p meta -q >> "$log" 2>&1; then
    echo -e "${master}\tPRODIGAL_FAIL\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    [[ "$KEEP_NEG" -eq 0 ]] && rm -f "$ntfa" "$aafa" "$gbk" "$tbl" 2>/dev/null || true
    return 0
  fi

  if ! hmmsearch --cpu "$CPU" -E "$EVALUE" --tblout "$tbl" -o /dev/null "$PROFILE" "$aafa" >> "$log" 2>&1; then
    echo -e "${master}\tHMMSEARCH_FAIL\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    [[ "$KEEP_NEG" -eq 0 ]] && rm -f "$ntfa" "$aafa" "$gbk" "$tbl" 2>/dev/null || true
    return 0
  fi

  read -r hits bestE bestS < <(summarize_tbl "$tbl")

  if [[ "$hits" -eq 0 ]]; then
    echo -e "${master}\tNO_HITS\t0\tNA\tNA" > "$OUTDIR/parts/summary.${master}.tsv"
    [[ "$KEEP_NEG" -eq 0 ]] && rm -f "$ntfa" "$aafa" "$gbk" "$tbl" 2>/dev/null || true
  else
    # hits: source_id = contig_id
    tblout_to_hits_wgs "$tbl" > "$OUTDIR/parts/hits.${master}.tsv"
    echo -e "${master}\tOK\t${hits}\t${bestE}\t${bestS}" > "$OUTDIR/parts/summary.${master}.tsv"
  fi
}

run_mode_wgs() {
  [[ -f "$INPUT" ]] || { echo "[ERROR] Modo WGS requiere -i/--input"; exit 1; }

  local list="$TMPDIR/wgs_masters.list"
  : > "$list"
  while IFS= read -r line; do
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    local id; id="$(extract_master_id "$line")"
    [[ -n "$id" ]] && echo "$id" >> "$list"
  done < "$INPUT"

  run_queue "$JOBS_WGS" "$list" process_one_wgs_master
}

# ---------------- SRA workers ----------------
prepare_sra_list_from_txt() {
  local infile="$1"
  local outfile="$2"

  [[ -f "$infile" ]] || { echo "[ERROR] No existe --sra-list: $infile"; return 2; }

  # Acepta ficheros TXT con un ID por línea. También tolera espacios/tabuladores,
  # líneas vacías y comentarios que empiezan por #. Se usa solo la primera columna.
  awk '
    BEGIN{IGNORECASE=0}
    {
      gsub(//, "", $0)
      if($0 ~ /^[[:space:]]*$/) next
      if($0 ~ /^[[:space:]]*#/) next
      print $1
    }
  ' "$infile" | grep -E '^(SRR|ERR|DRR)[0-9]+' | sort -u > "$outfile"

  [[ -s "$outfile" ]] || { echo "[ERROR] --sra-list no contiene IDs SRR/ERR/DRR válidos: $infile"; return 2; }
}

get_srrs_from_bioproject() {
  local prj="$1" outlist="$2"
  local runinfo="$TMPDIR/runinfo_${prj}.csv"
  esearch -db sra -query "$prj" | efetch -format runinfo > "$runinfo"
  local col
  col="$(head -n 1 "$runinfo" | awk -F',' '{for(i=1;i<=NF;i++) if($i=="Run"){print i; exit}}')"
  [[ -n "$col" ]] || { echo "[ERROR] No encuentro columna Run en runinfo"; return 2; }

  tail -n +2 "$runinfo" \
    | awk -F',' -v c="$col" '{if($c~/^(SRR|ERR|DRR)/) print $c}' \
    | sort -u > "$outlist"
}

save_spades_log_if_present() {
  local srr="$1"
  local asmdir="$2"
  local dest="$OUTDIR/logs/${srr}.spades.log"

  # SPAdes suele crear spades.log dentro del directorio de ensamblaje.
  # Lo copiamos antes de borrar intermedios para conservar el diagnóstico nativo.
  if [[ -s "$asmdir/spades.log" ]]; then
    cp "$asmdir/spades.log" "$dest" 2>/dev/null || true
    log_error "$srr: SPAdes falló; log nativo conservado en $dest"
  else
    log_error "$srr: SPAdes falló; no se encontró $asmdir/spades.log. Ver log combinado del SRR."
  fi
}

cleanup_sra_intermediates() {
  local srr="$1"
  local srrdir="$2"
  local readsdir="$3"
  local asmdir="$4"
  local ntfa="$5"
  local aafa="$6"
  local gbk="$7"
  local tbl="$8"

  # Limpieza agresiva de intermedios tras fallo de SPAdes/ensamblaje.
  # Se conservan únicamente logs y ficheros summary/hits ya generados.
  rm -rf "$srrdir" "$readsdir" "$asmdir" 2>/dev/null || true
  rm -f "$ntfa" "$aafa" "$gbk" "$tbl" 2>/dev/null || true
  log_warn "$srr: intermedios SRA eliminados tras fallo de ensamblaje/SPAdes"
}

process_one_srr() {
  local srr="$1"
  local log="$OUTDIR/logs/${srr}.sra.log"
  : > "$log"
  log_info "$srr: inicio procesamiento SRA"

  local srrdir="$OUTDIR/tmp/${srr}.prefetch"
  local readsdir="$OUTDIR/reads/$srr"
  local asmdir="$OUTDIR/assembly/$srr"
  local ntfa="$OUTDIR/nt/${srr}.contigs.fna"
  local aafa="$OUTDIR/aa/${srr}.contigs.faa"
  local gbk="$OUTDIR/nt/${srr}.prodigal.gbk"
  local tbl="$OUTDIR/hmm/${srr}.tblout"
  mkdir -p "$srrdir" "$readsdir" "$asmdir"

  if ! prefetch "$srr" -O "$srrdir" >> "$log" 2>&1; then
    echo -e "${srr}	PREFETCH_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_error "$srr: fallo en prefetch. Ver $log"
    rm -rf "$srrdir" 2>/dev/null || true
    return 0
  fi

  if ! fasterq-dump --threads "$SRA_THREADS" --split-3 --outdir "$readsdir" "$srrdir/$srr" >> "$log" 2>&1; then
    echo -e "${srr}	FASTQ_DUMP_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_error "$srr: fallo en fasterq-dump. Ver $log"
    [[ "$KEEP_READS" -eq 0 ]] && rm -rf "$readsdir" || true
    rm -rf "$srrdir" 2>/dev/null || true
    return 0
  fi

  # Preferir .gz si existen (si tú los introduces)
  local r1 r2 r3 s1
  r1="$(pick_read_file "$readsdir/${srr}_1")"
  r2="$(pick_read_file "$readsdir/${srr}_2")"
  r3="$(pick_read_file "$readsdir/${srr}_3")"
  s1="$(pick_read_file "$readsdir/${srr}")"

  local args=()
  if [[ -n "$r1" && -n "$r2" ]]; then
    args+=( --pe1-1 "$r1" --pe1-2 "$r2" )
    [[ -n "$r3" ]] && args+=( --pe1-s "$r3" )
  elif [[ -n "$s1" ]]; then
    args+=( --s1 "$s1" )
  elif [[ -n "$r1" && -z "$r2" ]]; then
    args+=( --s1 "$r1" )
  else
    echo -e "${srr}	READS_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_error "$srr: no se encontraron reads válidas tras fasterq-dump. Ver $log"
    [[ "$KEEP_READS" -eq 0 ]] && rm -rf "$readsdir" || true
    rm -rf "$srrdir" 2>/dev/null || true
    return 0
  fi

  local spades_rc=0
  log_info "$srr: inicio SPAdes assembler=$ASSEMBLER threads=$SPADES_THREADS memory=${SPADES_MEM}GB"
  if [[ "$ASSEMBLER" == "metaspades" ]]; then
    metaspades.py -o "$asmdir" "${args[@]}" --threads "$SPADES_THREADS" --memory "$SPADES_MEM" >> "$log" 2>&1 || spades_rc=$?
  else
    spades.py --meta -o "$asmdir" "${args[@]}" --only-assembler --threads "$SPADES_THREADS" --memory "$SPADES_MEM" >> "$log" 2>&1 || spades_rc=$?
  fi

  if [[ "$spades_rc" -ne 0 ]]; then
    echo -e "${srr}	SPADES_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    save_spades_log_if_present "$srr" "$asmdir"
    log_error "$srr: SPAdes terminó con código $spades_rc. Log combinado: $log"
    cleanup_sra_intermediates "$srr" "$srrdir" "$readsdir" "$asmdir" "$ntfa" "$aafa" "$gbk" "$tbl"
    return 0
  fi

  local contigs="$asmdir/contigs.fasta"
  if [[ ! -s "$contigs" ]]; then
    echo -e "${srr}	ASSEMBLY_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    save_spades_log_if_present "$srr" "$asmdir"
    log_error "$srr: SPAdes terminó sin error, pero no se generó contigs.fasta válido. Ver $log"
    cleanup_sra_intermediates "$srr" "$srrdir" "$readsdir" "$asmdir" "$ntfa" "$aafa" "$gbk" "$tbl"
    return 0
  fi

  cp "$contigs" "$ntfa"

  if ! prodigal -i "$ntfa" -a "$aafa" -o "$gbk" -p meta -q >> "$log" 2>&1; then
    echo -e "${srr}	PRODIGAL_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_error "$srr: fallo en prodigal. Ver $log"
    return 0
  fi

  if ! hmmsearch --cpu "$CPU" -E "$EVALUE" --tblout "$tbl" -o /dev/null "$PROFILE" "$aafa" >> "$log" 2>&1; then
    echo -e "${srr}	HMMSEARCH_FAIL	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_error "$srr: fallo en hmmsearch. Ver $log"
    return 0
  fi

  read -r hits bestE bestS < <(summarize_tbl "$tbl")

  if [[ "$hits" -eq 0 ]]; then
    echo -e "${srr}	NO_HITS	0	NA	NA" > "$OUTDIR/parts/summary.${srr}.tsv"
    rm -f "$ntfa" "$aafa" "$gbk" "$tbl" 2>/dev/null || true
    [[ "$KEEP_ASM" -eq 0 ]] && rm -rf "$asmdir" || true
    log_info "$srr: sin hits; intermedios eliminados según opciones"
  else
    tblout_to_hits_sra "$tbl" "$srr" > "$OUTDIR/parts/hits.${srr}.tsv"
    echo -e "${srr}	OK	${hits}	${bestE}	${bestS}" > "$OUTDIR/parts/summary.${srr}.tsv"
    log_info "$srr: OK hits=$hits bestE=$bestE bestS=$bestS"
  fi

  [[ "$KEEP_READS" -eq 0 ]] && rm -rf "$readsdir" || true
  rm -rf "$srrdir" 2>/dev/null || true
  log_info "$srr: fin procesamiento SRA"
}

run_mode_sra() {
  if [[ -z "$BIOPROJECT" && -z "$SRA_LIST" ]]; then
    echo "[ERROR] Modo SRA requiere --bioproject o --sra-list"
    exit 1
  fi

  if [[ -n "$BIOPROJECT" && -n "$SRA_LIST" ]]; then
    echo "[ERROR] Usa solo una fuente SRA: --bioproject o --sra-list, no ambas"
    exit 1
  fi

  local list="$TMPDIR/sra_runs.txt"

  if [[ -n "$SRA_LIST" ]]; then
    log_info "Preparando lista SRA desde TXT: $SRA_LIST"
    prepare_sra_list_from_txt "$SRA_LIST" "$list"
  else
    log_info "Obteniendo SRR desde BioProject: $BIOPROJECT"
    get_srrs_from_bioproject "$BIOPROJECT" "$list"
  fi

  if [[ "$MAX_RUNS" -gt 0 ]]; then
    head -n "$MAX_RUNS" "$list" > "$list.sub"
    mv "$list.sub" "$list"
  fi

  local n_sra
  n_sra="$(wc -l < "$list" | tr -d ' ')"
  log_info "SRA runs a procesar: $n_sra"

  run_queue "$JOBS_SRA" "$list" process_one_srr
}

# ---------------- merge parts ----------------
merge_outputs() {
  echo -e "$HITS_HDR" > "$HITS"
  echo -e "$SUM_HDR" > "$SUMMARY"

  # hits (si no existen, no pasa nada)
  if compgen -G "$OUTDIR/parts/hits.*.tsv" > /dev/null; then
    cat "$OUTDIR"/parts/hits.*.tsv >> "$HITS"
  fi
  # summary
  if compgen -G "$OUTDIR/parts/summary.*.tsv" > /dev/null; then
    cat "$OUTDIR"/parts/summary.*.tsv >> "$SUMMARY"
  fi
}

# ---------------- dispatch ----------------
log_info "Inicio pipeline MODE=$MODE OUTDIR=$OUTDIR PROFILE=$PROFILE EVALUE=$EVALUE CPU=$CPU JOBS=$JOBS JOBS_WGS=$JOBS_WGS JOBS_SRA=$JOBS_SRA"

if [[ "$MODE" == "WGS" ]]; then
  log_info "Ejecutando modo WGS"
  run_mode_wgs
else
  log_info "Ejecutando modo SRA BIOPROJECT=$BIOPROJECT SRA_LIST=$SRA_LIST MAX_RUNS=$MAX_RUNS ASSEMBLER=$ASSEMBLER"
  run_mode_sra
fi

log_info "Fusionando salidas"
merge_outputs

log_info "Pipeline finalizado OK. Hits=$HITS Summary=$SUMMARY"
echo "[INFO] OK"
echo "[INFO] Hits: $HITS"
echo "[INFO] Summary: $SUMMARY"
echo "[INFO] Log general: $PIPELINE_LOG"
