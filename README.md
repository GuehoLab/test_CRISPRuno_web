# CRISPRuno Web

[CRISPRuno](https://github.com/clementlab/CRISPRuno) 0.1.22(단방향 앵커 증폭 시퀀싱 분석)를 브라우저에서 실행하는 웹 버전입니다.
서버가 필요 없고, FASTQ·reference·결과 파일은 모두 이 컴퓨터의 브라우저 저장소(OPFS)에 남습니다.

- **bowtie2 대신 minibwa**(WebAssembly; wasm32 / wasm64 / 멀티스레드 wasm64)로 정렬합니다.
- **Cas-OFFinder 검색은 빠져 있습니다**(settings의 `PAM`, `casoffinder_num_mismatches`는 읽기만 하고 무시).
- CRISPResso2는 CRISPRuno가 쓰는 **정량 부분만** 이식했습니다(CRISPResso 사이트별 HTML/그림 없음).
- 그림은 matplotlib PNG/PDF 대신 **SVG**로 다시 그리고, 원본과 같은 구성의 **HTML 리포트**(`<root>.html`, 그림 내장, 오프라인으로 열림)를 만듭니다.
- CRISPRunoBatch / CRISPRunoCompare는 포함하지 않습니다.

## 실행

```
python serve.py            # http://localhost:8000  (COOP/COEP 헤더를 보내 멀티스레드 사용 가능)
```

1. **인덱스 만들기**(`build.html`): reference FASTA → `.l2b` + `.mbw`. 한 번 만들어 파일로 저장해 두면 계속 씁니다.
   CLI `minibwa index ref.fa prefix`로 만든 파일도 그대로 쓸 수 있습니다. (같은 origin에서 연 다른 minibwa 웹 앱과 OPFS의 인덱스를 공유합니다.)
2. **CRISPRuno**(`crispruno.html`)
   - settings: CRISPRuno settings 파일 형식(`key<TAB>value`) 그대로. 명령줄 옵션도 같은 형식으로 적습니다.
     파일 경로는 이름만 보고, 아래에서 고른 파일 중에서 찾습니다(`fastq_r1/fastq_r2/fastq_umi/cut_region_annotation_file`).
   - 입력: FASTQ(.gz 가능), minibwa 인덱스(브라우저 저장소 또는 `.l2b`+`.mbw` 파일), reference FASTA(+`.fai`; 없으면 `.l2b`에서 서열을 읽음).
   - 단계별 실행(체크한 단계만, 순서대로). 각 단계는 앞 단계가 run 폴더에 남긴 파일을 읽으므로 한 단계만 다시 돌릴 수 있습니다.

| 단계 | 원본 함수 (CRISPRunoCore.py) | 웹 버전 |
|---|---|---|
| prep | `prep_input` | primer·guide 서열을 minibwa로 정렬해 위치를 찾음(`-k 15 -m 15 -s 10`; bowtie2 `--end-to-end -c` 대신). origin 서열, read 수·평균 길이 |
| umi | `add_umi_from_umi_file`, `analyze_UMIs_initial` | 같음 (GINI는 정렬된 누적합으로 O(n log n), 값은 같음) |
| primer | `filter_on_primer` | SSW(`ssw.c`를 WebAssembly로 컴파일 — 같은 코드), cutadapt 3′ adapter 검색을 JS로 이식(k-mer 필터 + `Aligner.locate`) |
| align | `align_reads` | minibwa `-s 20 -m 20 --end-bonus=5 -b MD`(바꿀 수 있음; 아래 "정렬 옵션"). `samtools view -F 256 \| samtools sort`를 재현(보조 정렬 0x800 제거, `=` 표기, (tid,pos,strand) 안정 정렬) |
| assign | `make_final_read_assignments`, `collapse_cut_points` | 같음 (모든 `.final.*` 텍스트 파일) |
| crispresso | `prep_crispresso2`, CRISPResso2, `run_and_aggregate_crispresso` | CRISPResso2 2.3.4의 정량 경로 이식(아래) |
| summary | `make_final_summary`, `make_report` | 요약 + SVG 그림 + HTML 리포트 |

## 원본과의 동일성

같은 정렬(SAM)이 들어가면 나머지 단계의 출력 파일은 원본과 같도록 이식했습니다. `test/`의 회귀 테스트로 확인합니다:

- CRISPRuno 자체 테스트 18개 케이스(Cas-OFFinder 옵션 제외)를 Python(bowtie2)과 이 포트(같은 bowtie2 정렬)로 돌려 텍스트 출력 33–35개 파일이 모두 바이트 단위로 같음
  (`.final_assignments.txt`, `.final_cut_points*.txt`, `.cut_homology.txt`, `.cut_classification.txt`, 각종 histogram, `.CRISPResso.summary.txt`, `.annotated_final_assignments`, `.summary.txt` 등).
- SSW 트리밍 12,000 read, cutadapt `BackAdapter.match_to` 60,000 read, `global_align` + `find_indels_substitutions` 6,000 쌍: 차이 0.
- 실제 데이터(UDiTaS 1450-R-CL03-DE, 614,470 read; UMI, `dedup_by_final_cut_assignment_and_UMI`, `run_crispresso_on_novel_sites`, `write_discarded_read_info` 포함)를
  같은 bowtie2 정렬로 돌려 34/35 파일이 같음(다른 하나는 `commands.txt`의 절대 경로뿐).
- 브라우저 안의 WASM minibwa와 네이티브 minibwa(같은 옵션)의 결과: 21개 이상 파일이 같고 나머지는 read 1–2개 차이(같은 점수의 정렬 위치 선택 등).

정렬기만 다르므로 minibwa 결과는 bowtie2와 read 위치·soft clip·MAPQ가 다를 수 있습니다. 특히 MAPQ 척도가 다릅니다(bowtie2 최대 42, minibwa 최대 60;
`min_alignment_quality_score` 기본값 30은 그대로 씁니다).

### 정렬 옵션 (`--end-bonus`)

minibwa는 read 끝까지 정렬을 늘리는 보너스(end bonus, 기본 10)가 커서 read 끝의 mismatch 몇 개를 soft-clip하지 않고 그대로 정렬합니다.
CRISPRuno는 "양끝 일치 길이"(`arm_min_matched_start_bases`, MD 태그)로 poor alignment를 거르므로, 이 차이가 그대로 `Poor alignment` 수로 나타납니다.
그래서 이 앱의 minibwa에는 `--end-bonus=INT` 옵션을 추가했고(`engine/wasm/end-bonus.patch`), 기본값을 `--end-bonus=5`로 두었습니다.

실제 데이터(UDiTaS 1450-R, 61만 read)에서 bowtie2와 비교:

| | bowtie2 (원본) | minibwa 기본(eb 10) | **minibwa eb 5 (기본값)** | minibwa eb 0 |
|---|---:|---:|---:|---:|
| Poor alignment | 51,205 | 75,488 | **43,041** | 37,168 |
| 최종 분석 read | 74,250 | 73,281 | **75,220** | 75,607 |
| Linear | 35,779 | 35,080 | **36,453** | 36,658 |
| Chimera | 38,463 | 37,965 | **38,529** | 38,710 |
| CRISPResso 분석 read | 29,516 | – | **29,556** | – |

### 원본의 동작을 그대로 둔 것 (알아 둘 점)

- `additional_cut_site_file`은 0.1.22에서 읽기만 하고 `prep_input`에 넘기지 않아 효과가 없습니다(그대로).
- `short_indel_length_cutoff`는 명령줄 옵션에만 있고 실제로는 항상 50입니다(그대로).
- read 이름에 공백이 있으면(예: `@NB501311:...:17138 1:N:0:ATCG...`) 0.1.22는 정렬용 이름의 공백을 `_`로 바꾸지만 CRISPResso용 read를 찾을 때는 바꾸지 않아,
  CRISPResso 분석 read가 0개가 되고 모든 read가 `(not analyzed)`가 됩니다. 웹 버전도 기본은 같고, "웹 전용 옵션 → read ID 공백 보정"을 켜면 바로잡습니다.
- `collapse_cut_points`가 리스트를 순회하며 지우는 부분(다음 원소를 건너뜀), 원점 inversion/deletion 그림이 리포트에서 빠지는 오타 등도 출력이 같도록 그대로 재현했습니다
  (단, inversion/deletion 그림은 웹 리포트에는 넣었습니다).
- CRISPResso가 실패하는 조건(정렬 read 0개, quantification window가 제외 구간과 겹침, guide가 amplicon에 없음 등)은 같은 조건에서 `NA`로 기록합니다.

- bowtie2 2.5.x는 길이 0인 R2 mate가 있으면 "fewer reads in file specified with -2"로 멈춥니다. 비교용 참조 실행(`test/py/run_reference.py`, `test/run_case.js`)은 빈 mate를 `N`으로 바꿔 돌렸습니다. minibwa는 빈 mate도 그대로 받습니다.

### 의도적으로 다른 것

- 중간 FASTQ는 gzip 없이 저장합니다(`.fq`; 이름의 `.gz`만 다름). BAM 대신 정렬된 SAM(`<root>.genomeAlignment.sam`).
- `.info` 파일의 plot 객체 줄은 웹 리포트용 JSON입니다(원본은 PNG/PDF 경로).
- SAM 줄 끝의 개행을 떼고 태그를 읽습니다. 원본은 `MD:Z`가 줄의 마지막 태그면 개행 때문에 오른쪽 일치 길이를 0으로 읽는데, bowtie2는 MD를 마지막에 쓰지 않지만 minibwa는 쓸 수 있습니다.
- `samtools faidx` 대신 FASTA(+`.fai`)를 직접 읽습니다(시작 좌표 < 1은 1로).
- Cas-OFFinder 없음.

## 검증 기록

`test/run_case.js`(Node.js; 같은 `cu/*.js` 모듈을 파일 시스템 위에서 실행)와 `test/e2e.py`(headless Chromium으로 페이지 전체 실행):

```
# 1) 원본 결과 만들기 (CRISPRuno 0.1.22 + CRISPResso2 2.3.4 + cutadapt 5.2, bowtie2/samtools 필요)
python test/py/run_reference.py CRISPRuno/tests REF_DIR
# 2) 웹 포트로 같은 케이스 실행 + 비교 (bowtie2로 같은 정렬; --minibwa=prefix 이면 minibwa로 정렬)
for c in REF_DIR/*; do node test/run_case.js $c OUT/$(basename $c) --quiet; done
# 단위 비교
python test/py/gen_ssw_cases.py 1 300 > ssw.json && node test/test_ssw.js ssw.json
python test/py/gen_cutadapt_cases.py 1 4000 > ca.json && node test/test_cutadapt.js ca.json
python test/py/gen_align_cases.py 1 3000 > al.json && node test/test_align.js al.json
# 브라우저 전체 실행
python test/e2e.py --settings s.txt --files r1.fq.gz,r2.fq.gz,umi.fq.gz --index g.l2b,g.mbw --fasta g.fa --out out
```

## 속도

실제 데이터(614,470 read pair, UMI) 전체 실행: Python CRISPRuno(bowtie2, 1 core) 약 19분 → 이 앱(headless Chromium, 2코어 VM, 정렬 2 스레드) 약 90초
(prep 5 s, UMI 20 s, primer 18 s, 정렬 28 s, 할당 14 s, CRISPResso 4 s, 리포트 0.4 s).

## 빌드

- minibwa WebAssembly(`minibwa.wasm`, `minibwa64.wasm`, `minibwa64-mt.*`): `engine/build.sh`(wasm32), `engine/build64.sh`(wasm64, 멀티스레드).
  minibwa-web(GUIDE-seq 웹 버전)과 같은 패치(`minibwa-wasm.patch`)에 `end-bonus.patch`를 더해 빌드합니다. `wasi-shim.js`, `align-worker.js`, `worker.js`는 minibwa-web과 같습니다.
- `cu/ssw.wasm`: `wasm/ssw/build.sh` (zig 0.13 또는 wasi-sdk; SSE2 → WASM SIMD128 compat 헤더).

## 구성

```
serve.py                  로컬 서버 (COOP/COEP)
web/index.html            처음
web/build.html, build.js  인덱스 만들기 (minibwa index)
web/crispruno.html, .js   CRISPRuno 페이지
web/cu-worker.js          단계 실행 worker (OPFS)
web/align-worker.js       minibwa 정렬 worker
web/cu/settings.js        parse_settings
web/cu/pipeline.js        processCRISPRuno를 단계로 나눈 것 (prep_input, make_final_summary 포함)
web/cu/umi.js             UMI 추가·분석, GINI
web/cu/ssw.js + ssw.wasm  SSW (primer/origin 트리밍)
web/cu/cutadapt.js        cutadapt 3' adapter (BackAdapter, Aligner, k-mer 필터)
web/cu/primer.js          filter_on_primer
web/cu/sam.js             samtools view/sort 재현
web/cu/assign.js          make_final_read_assignments, collapse_cut_points
web/cu/cralign.js         CRISPResso2 global_align, find_indels_substitutions, Python round
web/cu/crispresso.js      prep_crispresso2, CRISPResso2 정량, run_and_aggregate_crispresso
web/cu/plots.js, report.js SVG 그림, HTML 리포트
web/lib/                  textio(FASTQ/FASTA/gzip 읽기), pyfmt(Python float 표기), fflate
wasm/ssw/                 ssw.wasm 빌드 (zig/wasi-sdk, SSE2 compat 헤더)
engine/                   minibwa WebAssembly 빌드 스크립트와 패치
test/                     회귀 테스트 (Node.js, Python 참조, Playwright)
```

## 라이선스

- CRISPRuno와 CRISPResso2는 MGH의 학술·비상업 라이선스(EULA)입니다. 이 웹 포트의 `cu/` 중 CRISPRuno/CRISPResso2에서 옮긴 코드는 그 조건을 따릅니다
  (학술 연구용, 재배포 금지 — 공개 배포 전에 라이선스 보유자와 확인하세요).
- minibwa: MIT. SSW(`ssw.c`): MIT. cutadapt: MIT. fflate: MIT. coi-serviceworker: MIT. EDNAFULL: NCBI.
